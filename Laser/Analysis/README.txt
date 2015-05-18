depending on your tree input you need to use
EEstudies_continued or EEstudies_VectorNTuples_continued.

The files without continued are here to create the TGraphs from the trees (VPT/PN change vs. time per crystal), so in principle one would have only needed one stream of files.
Each macro did one part of the study, e.g. first part
EEstudies_VectorNTuples_continued.C
then
EEstudies_VectorNTuples_continued_part2.C
etc.

The final plots were done without
EEstudies_VectorNTuples_everything_continued_daily.C
which should do everything in one go, but without detailed initial studies.
e.g. what crystals behave strange and therefore needed to be vetoed, initial studies showing what kind of fit we want to do.
(Of course except of making the tgraphs.)

I used the EEstudies_continued macros mainly for cross checks I did using the orange LED.