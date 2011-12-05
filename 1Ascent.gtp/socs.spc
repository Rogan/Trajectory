#--Settings for Parameter Estimation
#--Default settings for all phases
#gesop_sope_file=model_dir/src/gesop_daten_20Hz_ver_p1.txt
#gesop_sope_weights_file=model_dir/src/weights_p1.txt
#gesop_sope_states=1,3,4,6,9,11,12,13
#gesop_sope_controls=2
#
#--Phase specific settings
#gesop_sope_file_2=model/src/gesop_daten_20Hz_ver_p2.txt
#gesop_sope_weights_file_2=model/src/weights_p2.txt
#gesop_sope_states_2=1,3,4,6,9,11,12,13
#gesop_sope_controls_2=2
#
#
#--Setting for the output level
#-- A Mathematical Optimal Control Interpretation
#-- B Right Hand Side Sparsity Pattern
#-- C Transcription Summary
#-- D SOCS I/O Parameters
#-- E Variable Grid Map
#-- F Constraint Grid Map
#-- G Optimal Control Scale Information
#-- H Index Sets and Sparsity Summary
#-- I Variables
#-- J Objective Function
#-- K Constraints
#-- L Control Analysis Constraint Error
#-- M Relative DAE Error
#-- N Dynamic Status Display
#-- O Dynamic Constraint Adjoints
#-- P Dynamic Variable Adjoints
#-- Q Control Analysis Lagrangian Errors
#-- R Control Analysis Optimality Error
#-- S Grid Refinement Summary
#recommended:		SOCOUT=a1b1c1d1e0f0g1h1i4j0k0l1m1n0o0p0q0r0s1
#previous: 		SOCOUT=a1b1c1d1e4f4g1h1i4j4k4l0m0n1o0p0q0r0s1
#GESOP default: 	socout=a0b0c1d0e0f0g1h1i0j0k1l0m0n1o0p0q0r0s1
#
#
#--Settings to avoid saving large hold arrays in scratch file
INCORE = 1
#
#
#--Settings for selection of transcription method 
#ITSWCH=16
#MTSWCH=2
#NSSWCH=1
#
#
#--Maximum number of function evaluations (GESOP default is 150,000)
#hhsnlp:maxnfe=150000
#
#
#--Size of CSTAT which specifies the size of the real workspace for 
#--the transformation to b-spline (GESOP default is 200000) 
gesop_cs=800000
#
#--MXDATA is the maximum number of discrete data values per phase
#MXDATA=0
#
#--MXPARM is the maximum number of parameters per phase
#MXPARM= 30
#
#--MXPCON is the maximum number of user defined constraints per phase
#MXPCON= 60
#
#--MXSTAT is the maximum number of dynamic variables per phase
#MXSTAT= 12
#
#--MXTERM is the maximum number of terms per phase
#MXTERM= 12
#
#
#
#stskl: line 975 of gsocs.adb