function ff_dropboxInstall

eval('! cd ~ && wget -O - "https://www.dropbox.com/download?plat=lnx.x86_64" | tar xzf -')
eval('! ~/.dropbox-dist/dropboxd')

end