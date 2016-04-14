#!/bin/bash

# +----------------------------------------------------------------------+
# |                                                                      |
# |  Set up Mac OS X to store temporary files in RAM rather than on disk.|
# |  Updated for Yosemite, 16 GB of RAM                                  |
# |                                                                      |
# |  By James Newell <jrnewell@github>                                   |
# |                                                                      |
# |  Originally by Ricardo Gameiro <http://blogs.nullvision.com/?p=357>  |
# |  Changes by Daniel Jenkins                                           |
# |     <http://blogs.nullvision.com/?p=357#comment-1140>                |
# |  Additional changes by Benjamin Chu <http://intechnicolor.net>       |
# +----------------------------------------------------------------------+

# if you want a different location for startup script, change this variable
SCRIPT_DIR=/usr/local/sbin
SCRIPT=ram-temp-disks
if [ ! -d "$SCRIPT_DIR" ]; then
    echo "$SCRIPT_DIR does not exist.  Either create the directory or change SCRIPT_DIR in this script..."
    exit 1
fi
cd "$SCRIPT_DIR"

cat << "EOF" | tee "$SCRIPT" > /dev/null
#!/bin/sh

# Create a RAM disk with same perms as mount point
RAMDisk() {
    local mount_pt="$1"
    local ramdisk_size=$(($2*1024*1024/512))

    # create ram disk
    echo "Creating RamFS for $mount_pt"
    local ramdisk_dev="$(hdiutil attach -drivekey system-image=yes -nomount -noautoopen ram://${ramdisk_size} | sed -e 's/[[:space:]]*$//')"
    
    # seems to help with preventing issue where some disks did not mount
    sleep 1

    # success
    if [ $? -eq 0 ]; then

        # Create HFS on the RAM volume.
        sudo newfs_hfs -v "$3" $ramdisk_dev

        # Store permissions from old mount point.
        eval `/usr/bin/stat -s $mount_pt`

        # Mount the RAM disk to the target mount point.
        sudo mount -t hfs -o noatime -o union -o nobrowse "$ramdisk_dev" "$mount_pt"

        # Restore permissions like they were on old volume.
        sudo chown $st_uid:$st_gid "$mount_pt"
        sudo chmod $st_mode "$mount_pt"
    fi
}

# 512 MB for temp directory
RAMDisk "/private/tmp" 512 "RAM Temp"

# 128 MB for run directory
RAMDisk "/var/run" 128 "RAM Run"

EOF
sudo chmod u+x,g+x,o+x "$SCRIPT"

cd /Library/LaunchDaemons
cat << EOF | sudo tee ram-temp-disks.plist > /dev/null
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
    <dict>
        <key>Label</key>
        <string>RAM Temp Disks</string>

        <key>ProgramArguments</key>
        <array>
            <string>${SCRIPT_DIR}/${SCRIPT}</string>
        </array>

        <key>RunAtLoad</key>
        <true/>
</dict>
</plist>
EOF

exit 0
