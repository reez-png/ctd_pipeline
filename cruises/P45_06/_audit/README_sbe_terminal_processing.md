# Step 1 SBE terminal processing

Generated on: 2026-06-17T13:22:07

## Purpose

This folder contains the internal audit trail for official Sea-Bird SBE Data Processing runs controlled from Python.

Python is used to copy files, build commands, run official Sea-Bird modules, record logs and save summaries.

## Processing design

The notebook runs the linear CTD sequence first:

```text
DatCnv > AlignCTD > WildEdit > CellTM > Filter > LoopEdit > Derive
```

After Derive, the workflow branches into four Bin Average products from the same unbinned derived CNV:

```text
1 m full cast
1 m downcast only
1 m upcast only
1 s full cast for SUNA matching
```

## Why a no-space working folder is used

The Sea-Bird command line modules can fail when input, output, configuration or PSA paths contain spaces.

This version copies working inputs to:

```text
C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\_sbe_work
```

## Active raw input folder

```text
C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L0\CTD
```

## PSA setup folder

```text
C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\psa
```

## SBE executable folder

```text
C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32
```

## Internal output audit folder

```text
C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\_audit
```

## L1 CTD deliverable folder

```text
C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L1\CTD
```

## Module sequence and branches

```text
name,description,exe,psa,enabled,input_kind,output_folder,output_suffix,needs_config,extra_args,psa_candidates,input_folder,input_suffix
01_datcnv,Data Conversion from raw Sea-Bird hex to engineering/oceanographic units (L1 readable full cast),DatCnvW.exe,01_datcnv.psa,True,raw,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L1\CTD,.cnv,True,[],,,
02_alignctd,Align CTD sensor timing (reduces salinity/oxygen spiking),AlignCTDW.exe,02_alignctd.psa,True,previous,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\Align_CTD,_al.cnv,False,[],,,
03_wildedit,Despike obvious outliers (GATED: Sea-Bird manual recommends before LoopEdit),WildEditW.exe,03_wildedit.psa,False,previous,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\Align_CTD,_al_we.cnv,False,[],,,
04_celltm,"Correct conductivity cell thermal mass (GATED: 'Question for Drew', slide 6)",CellTMW.exe,04_celltm.psa,False,previous,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\Align_CTD,_al_ctm.cnv,False,[],,,
05_filter,Low-pass filter pressure before LoopEdit (GATED: not in v2 PSA set),FilterW.exe,05_filter.psa,False,previous,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\Align_CTD,_al_filt.cnv,False,[],,,
06_loopedit,Mark pressure reversals and slow movement loops (ship heave),LoopEditW.exe,06_loopedit.psa,True,previous,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\Loop_Edit,_loop.cnv,False,[],,,
07_derive,"Derive oceanographic variables (salinity, density, etc.) - unbinned hub",DeriveW.exe,07_derive.psa,True,previous,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\Derived_Parameter,_der.cnv,False,[],,,
08b_binavg_1m_down,"1 m depth-bin average, downcast only (science profiles, slide 6)",BinAvgW.exe,08b_binavg_1m_down.psa,True,derived_unbinned,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\CTD_1m_down,_1m_down.cnv,False,[],,,
08d_binavg_1s_full,"1 s time-bin average, full cast for SUNA matching (slide 7)",BinAvgW.exe,08d_binavg_1s_full.psa,True,derived_unbinned,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\CTD_1s,_1s.cnv,False,[],,,
09b_asciiout_1m_down,ASCII export of 1 m downcast CNV,ASCII_OutW.exe,09_asciiout_semicolon.psa,False,,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\CTD_1m_down\asc_semicolon,_1m_down_ascii.asc,False,[],"['09b_asciiout_1m_down.psa', '09_asciiout_semicolon.psa', '09_asciiout.psa']",C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\CTD_1m_down,_1m_down.cnv
09d_asciiout_1s_full,ASCII export of 1 s full cast CNV for SUNA matching,ASCII_OutW.exe,09_asciiout_semicolon.psa,False,,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\CTD_1s\asc_semicolon,_1s_ascii.asc,False,[],"['09d_asciiout_1s_full.psa', '09_asciiout_semicolon.psa', '09_asciiout.psa']",C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\CTD_1s,_1s.cnv

```

## Module validation

```text
step,name,enabled,exe,exe_found,exe_path,canonical_psa,selected_psa_found,selected_psa_path,output_folder,output_suffix,input_suffix,status
1,01_datcnv,True,DatCnvW.exe,True,C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32\DatCnvW.exe,01_datcnv.psa,True,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\psa\01_datcnv.psa,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L1\CTD,.cnv,,OK
2,02_alignctd,True,AlignCTDW.exe,True,C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32\AlignCTDW.exe,02_alignctd.psa,True,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\psa\02_alignctd.psa,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\Align_CTD,_al.cnv,,OK
3,03_wildedit,False,WildEditW.exe,,,03_wildedit.psa,,,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\Align_CTD,_al_we.cnv,,SKIPPED_DISABLED
4,04_celltm,False,CellTMW.exe,,,04_celltm.psa,,,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\Align_CTD,_al_ctm.cnv,,SKIPPED_DISABLED
5,05_filter,False,FilterW.exe,,,05_filter.psa,,,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\Align_CTD,_al_filt.cnv,,SKIPPED_DISABLED
6,06_loopedit,True,LoopEditW.exe,True,C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32\LoopEditW.exe,06_loopedit.psa,True,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\psa\06_loopedit.psa,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\Loop_Edit,_loop.cnv,,OK
7,07_derive,True,DeriveW.exe,True,C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32\DeriveW.exe,07_derive.psa,True,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\psa\07_derive.psa,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\Derived_Parameter,_der.cnv,,OK
8,08b_binavg_1m_down,True,BinAvgW.exe,True,C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32\BinAvgW.exe,08b_binavg_1m_down.psa,True,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\psa\08b_binavg_1m_down.psa,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\CTD_1m_down,_1m_down.cnv,,OK
9,08d_binavg_1s_full,True,BinAvgW.exe,True,C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32\BinAvgW.exe,08d_binavg_1s_full.psa,True,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\psa\08d_binavg_1s_full.psa,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\CTD_1s,_1s.cnv,,OK
10,09b_asciiout_1m_down,False,ASCII_OutW.exe,,,09_asciiout_semicolon.psa,,,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\CTD_1m_down\asc_semicolon,_1m_down_ascii.asc,_1m_down.cnv,SKIPPED_DISABLED
11,09d_asciiout_1s_full,False,ASCII_OutW.exe,,,09_asciiout_semicolon.psa,,,C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\cruises\P45_06\L2\CTD\CTD_1s\asc_semicolon,_1s_ascii.asc,_1s.cnv,SKIPPED_DISABLED

```
