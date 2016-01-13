#!/bin/bash

old=.gitignore
new=$old.$$


find . \( -name .git -o -name .sconf_temp -o -name '*.so' \) -prune \
     -o -type f -a -perm /111 -exec file {} \; |\
    grep ELF |\
    cut -d: -f1 |\
    perl -ne 'BEGIN {$f=shift;$t=`cat "$f"`;@a=split(/\n/,$t);foreach(@a){print"$_\n";last if(/^###END/);}}s#^./##;print if !$already{$_}++;' "$old" > "$new"

mv "$new" "$old"
