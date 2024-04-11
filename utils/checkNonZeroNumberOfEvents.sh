#!/bin/bash

for f in $(ls $1/*.root); do
  #echo $f
  root -l -b -e "TFile*_file0 = TFile::Open(\""$f"\");
                 if ( ((TTree*)_file0->Get(\"tout\"))->GetEntries() == 0 )
                   std::cout<<\"\\n---------\\n\"<<_file0->GetName()<<\"\\n---------\\n\"" -q
done
