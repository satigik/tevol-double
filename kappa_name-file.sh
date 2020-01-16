#!/bin/bash
## If you don't enter the third parameter, the
## script will execute the second instruction
if [ -n "$3" ]
then
    mv ILk1D.wt ILk1D_BrGc-k$1.wt
    mv ILk.wt ILk_BrGc-k$1.wt
    mv IL1D.wt IL1D_BrGc-k$1-NeNk$2.wt
    mv IL.wt IL_BrGc-k$1-NeNk$2.wt
else
    mv ILk1D.wt ILk1D-k$1.wt
    mv ILk.wt ILk-k$1.wt
    mv IL1D.wt IL1D-k$1-NeNk$2.wt
    mv IL.wt IL-k$1-NeNk$2.wt
    mv Fek1D.wt Fek1D-k$1.wt
    mv Fek.wt Fek-k$1.wt
    mv Fe1D.wt Fe1D-k$1-NeNk$2.wt
    mv Fe.wt Fe-k$1-NeNk$2.wt
fi
