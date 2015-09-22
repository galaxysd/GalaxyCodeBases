## After Install

````bash
mkdir ~/office2016fonts
cd ~/office2016fonts
sudo su

for APP in Word Excel PowerPoint OneNote Outlook; do echo [$APP];mkdir $APP;mv /Applications/Microsoft\ ${APP}.app/Contents/Resources/Fonts/* $APP/; done

for APP in Word Excel PowerPoint OneNote Outlook; do echo [$APP];cd $APP;cfv -rr -C -t sha1;mv $APP.sha1 ..;cd ..;done

````

## After Update

````bash
for APP in Word Excel PowerPoint OneNote Outlook; do echo [$APP];mkdir $APP;ls /Applications/Microsoft\ ${APP}.app/Contents/Resources/Fonts/;done
````
