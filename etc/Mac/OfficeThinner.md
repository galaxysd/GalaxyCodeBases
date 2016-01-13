# Records for OfficeThinner

## Issues

### A for effort, but dangerous #7

https://github.com/goodbest/OfficeThinner/issues/7

#### pbowden-msft commented 8 days ago

Hi there,

I give you props for ingenuity, but this is a pretty dangerous solution to recommend to folks out there. I own the installer for Office 2016 for Mac at Microsoft, and there are technologies that we're releasing soon (including delta updaters) that will cause major issues with your solution - including the inability for users to launch apps that have been 'thinned' using your technique.

Additionally, any type of change to our app bundles will break the code signature on the app, which presents a security concern. Unfortunately, your solution is absolutely not supportable nor recommended.
I understand your intent behind creating such a program as there is significant duplication of resources between each app; however, that is intentional based on the design-pattern mandated by Apple for sandboxed apps.

Please feel free to contact me if you need more information.

Paul Bowden  
Apple Platform Experience  
Microsoft Corporation  

#### goodbest commented 8 days ago

Thanks for your advice, I will warn the users about this.

I love Office, and I want to use Office 2016 smoothly on a 128GB Mac.  
That's the initial goal of this script.

By the way, It's really exciting to see delta updaters soon.

#### pbowden-msft commented 8 days ago

Thanks for your understanding. I certainly don't want to dampen your creativity - I just want to make sure folks out there don't end up in a situation where their apps can no longer boot after applying updates.

#### brianjking commented 8 days ago

@pbowden-msft Does Microsoft have any plans in place to thin down the bloat natively?

If I've ran the thinning script I have to re-install Office from my 365 login, correct?

Thanks!

#### pbowden-msft commented 8 days ago

Hi Brian,

In Office 2011, we did used a shared resource architecture which is why the apps were much smaller. With the 2016 release, we've changed to sandboxed apps, which is Apple's mandated architecture for the Mac App Store. Unfortunately, sandboxed apps must have all their resources self-contained within the app bundle. The suite installer already uses a technique similar to @goodbest's that allows us to remove the duplicates from the download pkg, but because of code signature constraints, the installer duplicates the resources at install time.

We are actively working on some changes that will reduce the overall size of the packages, including download on demand fonts and proofing tools. For example, many East Asian fonts are tens of megabytes in size, and for customers who never render those glyphs, they can save space on disk by not having those fonts present. Same goes for dictionaries and grammar checkers. These changes will allow for our initial download pkg to be smaller, and resources can be downloaded on-demand when you use/need them. For things you never use, they won't get downloaded, so you'll get a disk space saving.

Hope this helps.

Paul.

#### pbowden-msft commented 8 days ago

Oh, @brianjking yes, if you wish to re-install, just download the latest Office 2016 suite from http://macadmins.software or the O365 Portal.

#### brianjking commented 7 days ago

@pbowden-msft Based on your comment above it seems like I should download the individual apps I want to use from http://macadmins.software/ instead of using the full suite installer from macadmins.software or the O365 poral site.

Is this correct? If so, do I first have to uninstall my existing versions?

Thanks!

#### pbowden-msft commented 7 days ago

@brianjking no, you can either download the suite or individual apps depending on your needs. Nope, you don't need to delete or uninstall - the installer will blindly overwrite whatever you have already.
