# Fills file md5 for NCBI SRA uploading page

set myVar to the clipboard
# display dialog myVar

set AppleScript's text item delimiters to " *"
tell application "Google Chrome" to activate
delay 1
repeat with p in paragraphs of myVar
	set ns to text items of p
	#display dialog item 1 of ns
	tell application "System Events"
		keystroke item 2 of ns
		keystroke tab
		keystroke item 1 of ns
		keystroke tab
		keystroke tab
	end tell
	delay 3.5
end repeat
#display dialog "Done !"

#tell application "Google Chrome" to activate
#delay 1
#tell application "System Events"
#	keystroke "lane1_Undetermined_L001_210210.A2.GZXJ008.1.addN.gz.101.fq.gz"
#	keystroke tab
#	keystroke "e0ddffd4dc74f5309819400e24cfacf7"
#	keystroke tab
#	keystroke tab
#end tell
