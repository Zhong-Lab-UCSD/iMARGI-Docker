# Guides for Issues on Windows and macOS System

iMARGI-Docker can work on any Docker supported system. Generally, We recommend to use it on Linux system, because
Docker supports Linux system very well. If you are using it on Windows or macOS, you might encounter some small Docker
configuration problems. As you might not familiar with settings of Docker Desktop or Docker Toolbox, so we collected some
potential issues you might encountered and guide you to solve them.

- [Guides for Issues on Windows and macOS System](#guides-for-issues-on-windows-and-macos-system)
  - [Checklist of Docker Settings](#checklist-of-docker-settings)
  - [Issues](#issues)
    - [Cannot Install Docker Desktop](#cannot-install-docker-desktop)
    - [Cannot find files](#cannot-find-files)
    - [Connection issue](#connection-issue)
    - [KeyError: 'chr1'](#keyerror-chr1)

## Checklist of Docker Settings

**If you encounter any problem, please check these settings first, they can fix all the known problems.**

- Install the latest stable version of Docker
  
  First of all, you need to check your system. If your system matches the requirements of Docker Desktop, please use
  Docker Desktop, which is much better than Docker Toolbox. The latest Docker Desktop can be downloaded from Docker
  official web page. However, the latest Docker Toolbox needs to be downloaded from its GitHub repository.
  
  - Docker Desktop for Windows [Download Link](https://hub.docker.com/editions/community/docker-ce-desktop-windows?tab=description)
  
  - Docker Desktop for macOS [Download Link](https://hub.docker.com/editions/community/docker-ce-desktop-mac)

  - Docker Toolbox GitHub page [Download Link](https://github.com/docker/toolbox/releases)

- Enable Virtualization
  
  Docker needs CPU virtualization technology to support Windows and macOS. Most of modern computers support
  virtualization. For Apple computer, it's enabled as default. For Windows computer, you might need to enable it in BIOS.
  [Check here to see how to enable it in BIOS.](https://www.isumsoft.com/computer/enable-virtualization-technology-vt-x-in-bios-or-uefi.html)
  
  Besides, you also need to turn on Hyper-V if you want to use Docker Desktop for Windows 10 Pro, Enterprise or
  Education version. 
  [Check here to see how to turn on Hyper-V.](https://docs.microsoft.com/en-us/virtualization/hyper-v-on-windows/quick-start/enable-hyper-v)

- Start Docker
  
  After installing Docker, you need to start it as a backend engine. You can find a Docker icon (whale) in the task bar.
  You can learn the basic usage of Docker Desktop and Docker Toolbox from official docs.

  - Instructions of Docker Desktop for Windows [Instructions](https://docs.docker.com/docker-for-windows/)
  
  - Instructions of Docker Desktop for macOS [Instructions](https://docs.docker.com/docker-for-mac/)
  
  - Instructions of Docker Toolbox for Windows [Instructions](https://docs.docker.com/toolbox/toolbox_install_windows/)
  
  - Instructions of Docker Toolbox for macOS [Instructions](https://docs.docker.com/toolbox/toolbox_install_mac/)

- Memory Limit Setting
  
  Enough memory is important to the iMARGI data processing pipeline. The required amount of memory depends on the size
  of reference genome. For human genome, at least 8GB free memory are required by BWA. Hence,** the memory on the machine
  needs to be more than 8 GB, which usually is 16 GB**. If the memory is not enough, BWA will generate an empty BAM file,
  then it will throw out some strange error information in following steps, such as `"KeyError: 'chr1'"`.
  
  - **System memory requirement: 16 GB.**
  
  - If you are using Docker Desktop for Windows or macOS, you can easily change the settings by right click the
  Docker icon (Whale) in the task bar, then go to _Settings_ -> _Advanced_ to change memory and CPU limits.

  - If you are using Docker Toolbox for Windows or macOS, which uses VirtualBox as backend, so you need to open _VirtualBox_,
  then stop _default VM_, select it and click on _settings_, then make changes as you want.

- Shared Folder
  
  When we use iMARGI-Docker, we need to 'mount' the data analysis working directory to the container using `-v` option.
  Usually, it will be automatically mounted when you use `-v` option. However, if you encounter some problems similar
  to "file not found", you need to check the shared folder in the Docker Desktop setting panel or VirtualBox settings.
  Besides, in Windows, the path of folder is different to Linux. For example, `C:/Users/test space/imargi_example` needs
  to be written as `-v "/c/Users/test space/imargi_example":/imargi`.

## Issues

### Cannot Install Docker Desktop

- CPU virtualization enabled?
- Hyper-V enabled for Docker Desktop? Only Windows 10 Pro, Enterprise or Education version.
- Windows 10 Home version? Pleae use Docker Toolbox.

[Check detail information in Technical Notes](https://sysbio.ucsd.edu/imargi_pipeline/technical_note.html#install-docker-on-different-systems)

### Cannot find files

- Is that folder or drive shared in Docker? Check Docker Desktop setting panel or VirtualBox settings.
- Is the path correct? For example, `C:/Users/test space/imargi_example` needs to be written as
  `-v "/c/Users/test space/imargi_example":/imargi`.

### Connection issue

This issue was reported in an old version of Docker Toolbox on Windows. The error information might be:

```
wsarecv: A connection attempt failed because the connected party did not properly respond after a period of time, or established connection failed because connected host has failed to respond.
```

We found that it might be a bug of old Docker Toolbox. So update to the newest version will fix the problem, such as
[version 18.09.3](https://github.com/docker/toolbox/releases).

### KeyError: 'chr1'

It's caused by insufficient memory and BWA. The memory required by BWA depends on the size of the reference genome. For
the human genome, it needs about 6 GB free memory for building or loading index files. If the memory is not enough,
BWA will generate an empty BAM file without any warning or error information. Then in the following steps of the
pipeline will throw a lot of error information containing `KeyError 'chr1'`.

So, generally, the system should have more than 8 GB total physical memory, usually it's 16 GB.

Besides, there is memory limit setting of Docker on Windows and Mac. The default memory limit is 2 GB. You need to
increase it to 8 GB. [Check detail information in Technical Notes](https://sysbio.ucsd.edu/imargi_pipeline/technical_note.html#change-docker-memory-settings-on-windows-and-macos)
