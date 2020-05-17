---
title: Sherpa run-time settings
date: Sun May 17 20:21:55 2020
...

Note that parameters that are never accessed by Sherpa during its run will not be listed below. On the other hand, "accessed" does not necessarily mean that the parameter had any effect on the run.

In rare cases, an alternative default value is being used. These alternatives will be separated by "`-- AND --`" from the standard default, which will always be listed on top.

Customised settings
-------------------
Note that parameters that can take on different values because they are set within a list, for example `param: [{x: 1}, {x: 2}, ...]`, will not appear in the config-file or command-line columns. They will be listed in the final-value column, with each different value separated by an "`-- AND --`" line.

| parameter | default value | override by SHERPA | Sherpa.yaml | command line | final value |
|-|-|-|-|-|-|
ANALYSIS|  |  | Rivet |  | Rivet |  |
ANALYSIS\_OUTPUT| /home/hiro/Documents/Projects/UNI/Bachelor/prog/runcards/pp\_phaeno/runcards/basic/Analysis/ |  |  | analysis | analysis |  |
BEAMS| 0, 0 |  | 2212, 2212 |  | 2212, 2212 |  |
BEAM\_ENERGIES| 0, 0 |  | 6500 |  | 6500 |  |
BEAM\_REMNANTS| 1 |  | false |  | 0 |  |
EVENTS| 100 |  | 10000000/8 | 1000 | 1000 |  |
EVT\_OUTPUT| 3 |  |  | 0 | 0 |  |
EW\_SCHEME| alphamZ |  | alpha0 |  | alpha0 |  |
FRAGMENTATION| Ahadic |  | None |  | None |  |
ME\_GENERATORS| Comix, Amegic, Internal |  | Comix, Amegic |  | Comix, Amegic |  |
MI\_HANDLER| Amisic |  | None |  | None |  |
OUTPUT| 2 |  | 3 |  | 3 |  |
PDF\_LIBRARY|  |  | LHAPDFSherpa |  | LHAPDFSherpa |  |
PDF\_SET|  |  | NNPDF31\_lo\_as\_0118 |  | NNPDF31\_lo\_as\_0118 |  |
PDF\_SET\_VERSIONS|  |  | 0, 0 |  | 0, 0 |  |
PROCESSES:94 \-94 \-\> 22 22:Integration\_Error|  |  |  |  | 0\.001 |  |
PROCESSES:94 \-94 \-\> 22 22:Order:EW| \-1 |  |  |  | 2 |  |
PROCESSES:94 \-94 \-\> 22 22:Order:QCD| \-1 |  |  |  | 0 |  |
PSI:ITMIN| \-\- AND \-\-<br />1000 |  |  |  | 1000 |  |
PSI:MAXOPT| \-\- AND \-\-<br />3 |  |  |  | 3 |  |
PSI:NOPT| \-\- AND \-\-<br />7 |  |  |  | 7 |  |
RESULT\_DIRECTORY| /home/hiro/Documents/Projects/UNI/Bachelor/prog/runcards/pp\_phaeno/runcards/basic/Results/ |  | out/Results | Results | Results |  |
RIVET:ANALYSES|  |  | MC\_DIPHOTON\_PROTON |  | MC\_DIPHOTON\_PROTON |  |
RIVET:IGNOREBEAMS| 0 |  | 1 |  | 1 |  |
RIVET:USE\_HEPMC\_NAMED\_WEIGHTS| 1 | 0 |  |  | 0 |  |
SCALES| METS\{MU\_F2\}\{MU\_R2\}\{MU\_Q2\} |  | VAR\{1/2\*Abs2\(p\[0\]\+p\[1\]\)\} |  | VAR\{1/2\*Abs2\(p\[0\]\+p\[1\]\)\} |  |
SELECTORS|  |  | Eta, 22, \-2\.5, 2\.5<br />PT, 22, 20, 200000000 |  | Eta, 22, \-2\.5, 2\.5<br />\-\- AND \-\-<br />PT, 22, 20, 200000000 |  |
SHOWER\_GENERATOR| CSS |  | None |  | None |  |
Settings kept at their default value
-------------------
| parameter | default value |
|-|-|
1/ALPHAQED\(0\)| 137\.03599976 |  |
ABS\_ERROR| 0 |  |
ALPHAQED\_DEFAULT\_SCALE| 0 |  |
ALPHAS:FREEZE\_VALUE| 1 |  |
ALPHAS:PDF\_SET\_VERSION| 0 |  |
ALPHAS:USE\_PDF| 0 |  |
ALPHAS\(MZ\)| 0\.118 |  |
AMEGIC:DEFAULT\_GAUGE| 1 |  |
AMEGIC:PARTIAL\_COMMIT| 0 |  |
ANALYSIS\_WRITEOUT\_INTERVAL| 1\.84467440737e\+19 |  |
AS\_FORM| Smooth |  |
BATCH\_MODE| 1 |  |
BEAM\_1| 0 |  |
BEAM\_2| 0 |  |
BEAM\_ENERGY\_1| 0 |  |
BEAM\_ENERGY\_2| 0 |  |
BEAM\_POLARIZATIONS| 0 |  |
BEAM\_SMAX| 1 |  |
BEAM\_SMIN| 1e\-10 |  |
BEAM\_SPECTRA| Monochromatic, Monochromatic |  |
BEAM\_SPECTRUM\_1| Monochromatic |  |
BEAM\_SPECTRUM\_2| Monochromatic |  |
BRH\_VMODE| 0 |  |
BUNCHES|  |  |
CHECK\_LIBLOCK| 0 |  |
CHECK\_POLES| 0 |  |
CHECK\_WEIGHT| 0 |  |
CITATION\_DEPTH| 1 |  |
CI\_OMODE| 1 |  |
CKM:Order| 0 |  |
CKM:Output| 0 |  |
COLOUR\_RECONNECTIONS| 1 |  |
COLOUR\_RECONNECTIONS:PMODE| 0 |  |
COLOUR\_RECONNECTIONS:Q\_0| 1 |  |
COLOUR\_RECONNECTIONS:RESHUFFLE| 0\.333333333333 |  |
COLOUR\_RECONNECTIONS:RESTRING| 0\.333333333333 |  |
COLOUR\_RECONNECTIONS:R\_0| 1 |  |
COLOUR\_RECONNECTIONS:etaQ| 0\.16 |  |
COLOUR\_RECONNECTIONS:etaR| 0\.16 |  |
COLOUR\_SCHEME| 1 |  |
COMIX:AEXP| 0\.9 |  |
COMIX:BMODE| 1 |  |
COMIX:ECMODE| 2 |  |
COMIX:ITMAX| 1000000 |  |
COMIX:ITMIN| 1000 |  |
COMIX:MFAC| 1 |  |
COMIX:N\_GPL| 3 |  |
COMIX:OMODE| 3 |  |
COMIX:PARTIAL\_COMMIT| 0 |  |
COMIX:PG\_MODE| 0 |  |
COMIX:PMODE| D |  |
COMIX:PS\_CHTH| 0\.01 |  |
COMIX:SEXP| 0\.75 |  |
COMIX:SRBASE| 1\.05 |  |
COMIX:STEXP| 0\.001 |  |
COMIX:TEXP| 0\.9 |  |
COMIX:THEXP| 1\.5 |  |
COMIX:TMODE| 1 |  |
COMIX:VINTS| 8 |  |
COMIX:VL\_MODE| 0 |  |
COMIX:VMODE| 1 |  |
COMIX:VSOPT| 1 |  |
COMIX:WF\_MODE| 0 |  |
COMIX:ZMODE| 0 |  |
COMIX\_DEFAULT\_GAUGE| 1 |  |
COMPRESS\_PARTONIC\_DECAYS| 1 |  |
CORE\_SCALE| Default |  |
COUPLINGS| Alpha\_QCD 1 |  |
DEBUG\_INTERVAL| 0 |  |
DEBUG\_STEP| \-1 |  |
DECAYER| 0 |  |
DECOMPOSE\_4G\_VERTEX| 1 |  |
DIPOLES:ALPHA| 1 |  |
DIPOLES:ALPHA\_FF| 1 |  |
DIPOLES:ALPHA\_FI| 1 |  |
DIPOLES:ALPHA\_IF| 1 |  |
DIPOLES:ALPHA\_II| 1 |  |
DIPOLES:AMIN| 1e\-08 |  |
DIPOLES:KAPPA| 0\.666666666667 |  |
DIPOLES:KT2MAX| 169000000 |  |
DIPOLES:NF\_GSPLIT| 5 |  |
ENHANCE\_XS| 0 |  |
ERROR| 0\.01 |  |
EVENT\_DISPLAY\_INTERVAL| 100 |  |
EVENT\_GENERATION\_MODE| PartiallyUnweighted |  |
EVENT\_INPUT|  |  |
EVENT\_OUTPUT|  |  |
EVENT\_SEED\_FILE| ran\.stat\.9872 |  |
EVENT\_SEED\_MODE| 0 |  |
EVENT\_TYPE| StandardPerturbative |  |
EVT\_FILE\_PATH| \. |  |
EVT\_OUTPUT\_START| 1 |  |
EW\_REN\_SCHEME| alpha0 |  |
EXTERNAL\_RNG| None |  |
FACTORIZATION\_SCALE\_FACTOR| 1 |  |
FINISH\_OPTIMIZATION| 1 |  |
FLAG\_PARTONIC\_DECAYS| 1 |  |
FREEZE\_PDF\_FOR\_LOW\_Q| 0 |  |
GENERATE\_RESULT\_DIRECTORY| 1 |  |
GLOBAL\_KFAC| 0 |  |
HARD\_DECAYS:Enabled| 0 |  |
HARD\_SPIN\_CORRELATIONS| 0 |  |
HELICITY\_SCHEME| 1 |  |
HEPMC\_EXTENDED\_WEIGHTS| 0 |  |
HEPMC\_INCLUDE\_ME\_ONLY\_VARIATIONS| 0 |  |
HEPMC\_TREE\_LIKE| 0 |  |
HEPMC\_USE\_NAMED\_WEIGHTS| 0 |  |
HISTOGRAM\_OUTPUT\_PRECISION| 6 |  |
IB\_THRESHOLD\_KILL| \-1e\+12 |  |
IB\_WHBINS| 100 |  |
INIT\_ONLY| 0 |  |
INTEGRATION\_ERROR| 0\.01 |  |
INTEGRATOR| Default |  |
INTRINSIC\_KPERP:CUT\_EXPO| 5 |  |
INTRINSIC\_KPERP:FORM| gauss\_limited |  |
INTRINSIC\_KPERP:MAX| 3 |  |
INTRINSIC\_KPERP:MEAN| 0 |  |
INTRINSIC\_KPERP:Q2| 0\.77 |  |
INTRINSIC\_KPERP:REFE| 7000 |  |
INTRINSIC\_KPERP:SCALE\_EXPO| 0\.08 |  |
INTRINSIC\_KPERP:SIGMA| 1\.5 |  |
INT\_MINSIJ\_FACTOR| 1e\-12 |  |
ISR\_SMAX| 1 |  |
ISR\_SMIN| 1e\-10 |  |
JET\_MASS\_THRESHOLD| 10 |  |
KFACTOR| None |  |
LHAPDF:DISALLOW\_FLAVOUR|  |  |
LHAPDF:USE\_Q2LIMIT| 1 |  |
LHEF\_PDF\_NUMBER| \-1 |  |
LOG\_FILE|  |  |
MASSIVE\_PS|  |  |
MASSLESS\_PS|  |  |
MAX\_TRIALS| 1000000 |  |
MCNLO\_DADS| 1 |  |
MEH\_EWADDMODE| 0 |  |
MEH\_NLOADD| 1 |  |
MEH\_QCDADDMODE| 0 |  |
MEMLEAK\_WARNING\_THRESHOLD| 16777216 |  |
MENLOPS\_MAX\_KFAC| 10 |  |
MEPSNLO\_PDFCT| 1 |  |
METS:CLUSTER\_MODE| 0 |  |
METS\_BBAR\_MODE| 1 |  |
ME\_QED:CLUSTERING\_ENABLED| 1 |  |
ME\_QED:CLUSTERING\_THRESHOLD| 10 |  |
ME\_QED:ENABLED| 1 |  |
ME\_QED:INCLUDE\_RESONANCES| 0 |  |
MI\_PDF\_LIBRARY|  |  |
MODEL| SM |  |
MPI\_EVENT\_MODE| 0 |  |
MPI\_OUTPUT| 0 |  |
MPI\_PDF\_SET|  |  |
MPI\_PDF\_SET\_VERSIONS|  |  |
MPI\_PT\_MAX| 1e\+12 |  |
MPI\_PT\_Max\_Fac| 1 |  |
MPI\_SEED\_MODE| 0 |  |
NLO\_IMODE| IKP |  |
NLO\_NF\_CONVERSION\_TERMS| None |  |
NLO\_SMEAR\_POWER| 0\.5 |  |
NLO\_SMEAR\_THRESHOLD| 0 |  |
NLO\_SUBTRACTION\_SCHEME| 0 |  |
NO\_ZERO\_PDF| 0 |  |
NUM\_ACCURACY| 1e\-10 |  |
ORDER\_ALPHAS| 2 |  |
OVERRIDE\_PDF\_INFO| 0 |  |
OVERWEIGHT\_THRESHOLD| 1e\+12 |  |
PB\_USE\_FMM| 0 |  |
PRETTY\_PRINT| On |  |
PRINT\_PS\_POINTS| 0 |  |
PRINT\_VERSION\_INFO| 0 |  |
PROCESSES:94 \-94 \-\> 22 22:CKKW|  |  |
PROCESSES:94 \-94 \-\> 22 22:Cut\_Core| 0 |  |
PROCESSES:94 \-94 \-\> 22 22:Decay|  |  |
PROCESSES:94 \-94 \-\> 22 22:DecayOS|  |  |
PROCESSES:94 \-94 \-\> 22 22:No\_Decay|  |  |
PSI:NPOWER| 1 |  |
PSI:NRAWMAX| 18446744073709551615 |  |
PSI:STOPOPT| 0 |  |
PSI:TIMESTEP\_OFFSET| 0 |  |
PSI:TIMESTEP\_SLOPE| 0 |  |
PS\_PT\_FILE|  |  |
Q2\_AS| 1 |  |
RANDOM\_SEED| \-1, \-1, \-1, \-1 |  |
RANDOM\_SEED1| \-1 |  |
RANDOM\_SEED2| \-1 |  |
RANDOM\_SEED3| \-1 |  |
RANDOM\_SEED4| \-1 |  |
REMNANTS:DELTA\_MASS| 1\.5 |  |
REMNANTS:SOFT\_ETA\_RANGE| 7\.5 |  |
REMNANTS:SOFT\_MASS| 5 |  |
REMNANTS:SOFT\_X\_EXPONENT| \-2 |  |
RENORMALIZATION\_SCALE\_FACTOR| 1 |  |
RESPECT\_MASSIVE\_FLAG| 0 |  |
REWEIGHT\_SPLITTING\_ALPHAS\_SCALES| 0 |  |
REWEIGHT\_SPLITTING\_PDF\_SCALES| 0 |  |
RIVET:\-l| 20 |  |
RIVET:HEPMC\_OUTPUT\_PRECISION| 15 |  |
RIVET:INCLUDE\_HEPMC\_ME\_ONLY\_VARIATIONS| 0 |  |
RIVET:JETCONTS| 0 |  |
RIVET:SKIPWEIGHTS| 1 |  |
RIVET:SPLITCOREPROCS| 0 |  |
RIVET:SPLITSH| 0 |  |
RIVET:USE\_HEPMC\_EXTENDED\_WEIGHTS| 0 |  |
RIVET:USE\_HEPMC\_SHORT| 0 |  |
RIVET:USE\_HEPMC\_TREE\_LIKE| 0 |  |
RIVET:XS\_OUTPUT\_PRECISION| 6 |  |
RLIMIT\_AS| 16615927808 |  |
RLIMIT\_BY\_CPU| 0 |  |
RUN\_MASS\_BELOW\_POLE| 0 |  |
SAVE\_STATUS|  |  |
SCALE\_FACTOR| 1 |  |
SELECTION\_WEIGHT\_MODE| 0 |  |
SHERPA\_CPP\_PATH|  |  |
SHERPA\_LDADD|  |  |
SHERPA\_LIB\_PATH|  |  |
SHERPA\_VERSION|  |  |
SHOW\_ANALYSIS\_SYNTAX| 0 |  |
SHOW\_FILTER\_SYNTAX| 0 |  |
SHOW\_ME\_GENERATORS| 0 |  |
SHOW\_MODEL\_SYNTAX| 0 |  |
SHOW\_NLOMC\_GENERATORS| 0 |  |
SHOW\_NTRIALS| 0 |  |
SHOW\_PDF\_SETS| 0 |  |
SHOW\_PS\_GENERATORS| 0 |  |
SHOW\_SCALE\_SYNTAX| 0 |  |
SHOW\_SELECTOR\_SYNTAX| 0 |  |
SHOW\_SHOWER\_GENERATORS| 0 |  |
SHOW\_VARIABLE\_SYNTAX| 0 |  |
SOFT\_COLLISIONS| None |  |
SOFT\_SPIN\_CORRELATIONS| 0 |  |
SP:ADD\_DOC| 0 |  |
SP:SET\_COLORS| 0 |  |
STATUS\_PATH|  |  |
THRESHOLD\_ALPHAS| 1 |  |
TIMEOUT| \-1 |  |
USERHOOKS| None |  |
USR\_WGT\_MODE| 1 |  |
VARIATIONS\_INCLUDE\_CV| 0 |  |
VEGAS\_MODE| 2 |  |
VIRTUAL\_EVALUATION\_FRACTION| 1 |  |
WIDTH\_SCHEME| CMS |  |
WRITE\_REFERENCES\_FILE| 1 |  |
YFS:1/ALPHAQED| 0 |  |
YFS:CHECK\_FIRST| 0 |  |
YFS:DRCUT| 1\.79769313486e\+308 |  |
YFS:FF\_RECOIL\_SCHEME| 2 |  |
YFS:FI\_RECOIL\_SCHEME| 2 |  |
YFS:INCREASE\_MAXIMUM\_WEIGHT| 1 |  |
YFS:IR\_CUTOFF| 0\.001 |  |
YFS:IR\_CUTOFF\_FRAME| Multipole\_CMS |  |
YFS:MAXEM| 2147483647 |  |
YFS:MINEM| 0 |  |
YFS:MODE| Full |  |
YFS:REDUCE\_MAXIMUM\_ENERGY| 1 |  |
YFS:STRICTNESS| 0 |  |
YFS:USE\_ME| 1 |  |
YFS:USE\_RUNNING\_PARAMETERS| 0 |  |
YFS:UV\_CUTOFF| 1\.79769313486e\+308 |  |
YUKAWA\_MASSES| Running |  |
