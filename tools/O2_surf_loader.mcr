#!MC 1120
#-----------------------------------------------------------
# Tecplot Macro for visualizing grids with associated fields
#-----------------------------------------------------------
#
# Get base name of files:
#
$!PromptforTextString |vbasename|
Instructions = "Enter the basename of the plot file [Space for default: O2_STEC]"
#
# 
#
$!IF "|vbasename|"==" "
$!VarSet |vbasename|="O2_STEC"
$!ENDIF
$!PromptforTextString |nvols|
Instructions = "Enter the number of files [Space for single file]"
$!VarSet |pvname|="|macrofilepath|/|vbasename|"
$!IF "|opt|"!="APPEND"
$!VarSet |opt|="NEW"
$!ENDIF
$!IF "|nvols|"==" "
$!READDATASET  '"|pvname|.plt"'
  READDATAOPTION = |opt|
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  INCLUDEDATA = YES
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  #VARPOSITIONLIST = <set> 
  #ZONELIST = <set> 
$!VarSet |opt|="APPEND"
$!ELSE
#------------------- Do loop for grid
$!LOOP |nvols|
$!READDATASET  '"|pvname||loop%d|.plt"'
  READDATAOPTION = |opt|
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  INCLUDEDATA = YES
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  #VARPOSITIONLIST = <set> 
  #ZONELIST = <set> 
$!VarSet |opt|="APPEND"
$!ENDLOOP
$!ENDIF
$!RemoveVar |nvols|
$!RemoveVar |pvname|
$!RemoveVar |vbasename|
