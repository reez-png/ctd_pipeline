# Processing SBE 55 Auto-Fire Logs (.afm) into a Bottle Table

**What this does:** turns the raw SBE 55 auto-fire log (`.afm`) for a cast into calibrated bottle
firing records (bottle number, time, pressure/depth, T, S at each firing). Those records are then
joined with the lab nitrate values to fill `<CRUISE>_bottle_nitrate.csv` for step03.

**Verified against:** Sea-Bird *SBE Data Processing* manual, rev. 7.26.7, Section 5 (pp. 73–82).

**The chain:**
```
.afm  +  .hex  --[Data Conversion]-->  .ros  (+ .bl auto-created)  --[Bottle Summary]-->  .btl
```
The `.afm` is not processed by a separate tool. It is a **scan-range source** for Data Conversion \u2014 the
same module (`DatCnvW.exe`) step01 already uses. Bottle Summary then reads the resulting `.ros` and
writes the human-readable `.btl` table.

---

## Prerequisites (per cast)

- The cast's raw CTD file: `L0\CTD\<CRUISE>_CTD_NN.hex`
- The cast's auto-fire log: `metadata\bottle_fire\<CRUISE>_CTD_NN.afm`
- The instrument config: `Instrument_Calibration_Files\19-8460.xmlcon`

**Important:** the `.afm` and the `.hex` must have the **same base file name** and be in the **same
directory** for Data Conversion to auto-match them. See the "File naming" note at the end \u2014 you will
likely need to co-locate a copy of the `.afm` next to the `.hex`, both named `<CRUISE>_CTD_NN`.

---

## Step 1 \u2014 Data Conversion (creates the .ros and .bl)

This is a **second, separate** Data Conversion run from the one in step01 (that one makes the profile
`.cnv`; this one makes the bottle `.ros`). Make a dedicated PSA so the two never collide.

Open **SBE Data Processing \u2192 Run \u2192 Data Conversion.**

**File Setup tab:**
1. **Instrument configuration file** \u2192 browse to `19-8460.xmlcon`.
2. **Input directory / files** \u2192 the cast `.hex` (the `.afm` must sit beside it, same base name).
3. **Output directory** \u2192 `metadata\bottle_fire\` (or a `bottle_ros\` subfolder if you prefer).
4. **Name append** \u2192 leave blank.

**Data Setup tab \u2014 the settings that matter:**
5. **Process scans to end of file** \u2192 checked.
6. **Convert data from** \u2192 **Upcast and downcast** (bottles fire on the upcast).
7. **Create file types** \u2192 select **"Create bottle file (.ros)"** (or "Create converted data and bottle
   file" if you also want a .cnv from this run \u2014 not required, step01 already made the profile .cnv).
8. **Source of scan range data for bottle file** \u2192 select **"ECO .afm file"**
   (the option reads roughly: *"Define scans from bottle fire module or ECO .afm file"*).
   - With this selected, Data Conversion reads the `.afm`, and **auto-creates a `.bl` file** (same base
     name, `.bl` extension) alongside the `.ros`. The `.bl` holds bottle sequence #, position, date,
     time, and scan numbers \u2014 Bottle Summary needs it in Step 2.
9. **Select Output Variables** \u2192 at minimum include **Pressure, Temperature, Conductivity (or
   Salinity)** \u2014 Bottle Summary requires these. Add **Depth [salt water, m]**, **Time**, and any others
   you want per-bottle (oxygen, fluorescence).

**Scan range offset / duration** (how many scans around each firing go into the `.ros`):
10. These define the averaging window per bottle. Sea-Bird's worked example uses **offset \u22122 s,
    duration 5 s**. For a 19plus V2 (~4 scans/s), that captures the firing plus a short window.
    A reasonable starting point: **Scan range offset = 0**, **Scan range duration = 2 s** (tight window
    right at firing), or the \u22122 s / 5 s example if you want more context. Note the max is 1440 scans
    per bottle (Bottle Summary limit).

11. **Save As** \u2192 `psa\bottle_01_datcnv_ros.psa` (a dedicated bottle-conversion PSA).
12. **Start Process.** Output: `<CRUISE>_CTD_NN.ros` and `<CRUISE>_CTD_NN.bl` in the output folder.

---

## Step 2 \u2014 Bottle Summary (creates the .btl)

Open **SBE Data Processing \u2192 Run \u2192 Bottle Summary.**

**File Setup tab:**
1. **Instrument configuration file** \u2192 `19-8460.xmlcon`.
2. **Input** \u2192 the `.ros` file from Step 1. (Bottle Summary automatically uses the matching `.bl` in the
   same folder for bottle position/time \u2014 keep them together.)
3. **Output directory** \u2192 `metadata\bottle_fire\`.
4. **Name append** \u2192 blank.

**Data Setup tab:**
5. **Averaged variables** \u2192 select **Pressure, Temperature, Salinity** (mean + std dev output per
   bottle). Add **Depth** if available.
6. **Derived variables** \u2192 optional (density, oxygen if you have it in the .ros).
7. **Output min/max values** \u2192 optional; not needed for the bottle table.

8. **Save As** \u2192 `psa\bottle_02_summary.psa`.
9. **Start Process.** Output: `<CRUISE>_CTD_NN.btl` \u2014 a plain-text table with, per bottle:
   **bottle position, date/time, and mean pressure / temperature / salinity (and depth).**

---

## Step 3 \u2014 From .btl to the bottle table (Python, later)

The `.btl` is plain text and easy to read. A small Python step (to be built and tested once you have a
real `.btl`) will:
1. Read each cast's `.btl` \u2192 bottle #, firing time, pressure\u2192depth, T, S.
2. Join the **lab nitrate** values by **bottle number** (per cast).
3. Write `metadata\bottle_nitrate\<CRUISE>_bottle_nitrate.csv` in the step03 schema
   (`cruise_id, cast_id, station, depth_m, time_utc, nitrate_uM, replicate_id, flag`).

step03 then consumes that CSV. **step03 never reads the `.afm`, `.ros`, `.bl`, or `.btl` directly** \u2014
the conversion is entirely a pre-step.

---

## File naming note (the one gotcha)

Data Conversion auto-matches the `.afm` to the `.hex` **only if they share a base name and folder**.
Your `.hex` is `<CRUISE>_CTD_NN.hex`, but the SBE 55 offload names the `.afm` however you saved it
(e.g. `afm.afm` on the sample we looked at). So before Step 1, for each cast:

- Copy the cast's `.afm` next to its `.hex`, renamed to match: `<CRUISE>_CTD_NN.afm` beside
  `<CRUISE>_CTD_NN.hex`.
- Keep the original in `metadata\bottle_fire\` as the archive; use the co-located copy for processing.

This is why the collection guide has you name the `.afm` `<CRUISE>_CTD_NN.afm` in the first place \u2014 it
makes this matching automatic.

---

## Could this be automated into step01?

Yes, eventually. Step 1 uses the same `DatCnvW.exe` step01 already drives, and Step 2's `BottleSumW.exe`
is scriptable the same way. Once the bottle PSAs exist and are tested on a real cast, the two runs
could be added as an optional bottle branch in step01's orchestration \u2014 so a single run produces both
the profile products and the `.btl` bottle summaries. Defer until tested on real `.afm` data.
