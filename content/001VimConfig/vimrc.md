# Vim 配置

```vimrc
syntax on
set nocompatible
set backspace=2
set tabstop=4
set softtabstop
set shiftwidth=4
set smarttab
set autoindent
set cindent
set smartindent
set number
set ruler
set cursorline
set scrolloff=3
set magic
set showmatch
set showcmd
set history=1000
set mouse=a
set hlsearch
set incsearch
set ignorecase
set smartcase

set printoption=paper:A4,syntax:n,wrap:y
set printdevice={{printer-name}}
noremap <C-p> :hardcopy<CR>
inoremap <C-p> <C-o>:hardcopy<CR>

map <C-a> gg0vG$

noremap <C-z> u
inoremap <C-z> <C-o>u

inoremap {<CR> {}<LEFT><CR><CR><UP><TAB>

map <C-p> i#include <bits/stdc++.h><CR><CR>int main() {<CR><CR>return 0;<UP><ESC>

map <F5> :w<CR>:!g++ % -o %< -lm -O2 -std=c++11; echo Exit code: $?; echo "==================="<CR>
imap <F5> <ESC>:w<CR>:!g++ % -o %< -lm -O2 -std=c++11; echo Exit code: $?; echo "==================="<CR>
```
