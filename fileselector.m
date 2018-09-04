% files = fileselector(dir,prefix,suffix)

function files = fileselector(dir,prefix,suffix)

files = spm_select('FPList', dir, sprintf('^%s.*%s',prefix,suffix));