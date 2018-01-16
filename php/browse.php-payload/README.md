# browse.php Payload

**browse.php** is a PHP single-file script providing an HTML based directory
browser. Once arbitrary write access is granted to a directory, deploying this
file will yield comprehensive insights in the directory structure of the server.

The script allows to...

* ... browse the file system with the privileges of the executing user
* ... download files using PHP's `readfile()`
* ... dictionary traversal, i.e. attempting paths, like `"..\..\..\"` *(if possible)*

## Demonstration

Let's have a look at some comprehensive screenshots. For demonstration purposes,
I deployed this script on my own server, which I obviously *have* write access
to. While taking the screenshots, I have picked a non-guessible name in order to
avoid "accidents" in the time beeing. After I was finished, I deleted the file -
**don't forget this if you test the script!**

### 1. Starting the script

You have done 99% of the work by acquiring write access remotely.
Congratulations, you are a genious! ;) Now, browse.php will list the current
directory. From there, navigation is a simple & UI based task.

![](https://bytecode77.com/images/sites/hacking/payloads/browse/001.png)

### 2. View & download files

You can view or download any file. Especially PHP files - and we all know which
ones are particularly interesting. This is a lot more convenient than the
`readfile("[...]\config.php")` code is that usually deployed in multiple trial
and error attempts until the correct path is hit.

![](https://bytecode77.com/images/sites/hacking/payloads/browse/002.png)

Any file you care about is accessible and can be downloaded. Note, that this is
a simple and therefore deployable & compatible script, not a feature complete
"remote cloud solution payload".

![](https://bytecode77.com/images/sites/hacking/payloads/browse/003.png)

Directory traversal attacks, like `"..\..\..\..\file.txt"` are possible, as we
can specify any path to the script. Here, I have deliberately weakened the
server configuration to demonstrate how a user that is not jailed could cause
you harm.

Please note, that I'm using **my own** server for this demonstration. I really
hope I didn't forget to delete the file afterwards...

![](https://bytecode77.com/images/sites/hacking/payloads/browse/004.png)

## Use cases

I have actually developed this while testing suphp, Apache MPM and similar, not
in an actual pentest. This helped me to debug through the web server
implementation. However, the main purpose primarily suits pentesting.

## Project Page

[![](https://bytecode77.com/images/shared/favicon16.png) bytecode77.com/hacking/payloads/browse](https://bytecode77.com/hacking/payloads/browse)