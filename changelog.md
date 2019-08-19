# Changelog
Visual Tracker (vistrack) change log
James Jun, Flatiron Institute

## [0.4.1] - 2019-08-19
### Fixed
- Bug fixed exporting `_posture.csv`
  - `_posture.csv` table contains following five feature points (from head to tail direction)
    '  Columns: x1(m), y1(m), x2(m), y2(m), x3(m), y3(m), x4(m), y4(m), x5(m), y5(m)',
    '    x1(m): x coordinate of the head tip @ grid frame of reference',
    '    y1(m): y coordinate of the head tip @ grid frame of reference',
    '    x2(m): x coordinate of the head-mid section @ grid frame of reference',
    '    x3(m): x coordinate of the mid section @ grid frame of reference',
    '    x4(m): x coordinate of the mid-tail section @ grid frame of reference',
    '    x5(m): x coordinate of the tail tip @ grid frame of reference',    

## [0.4.0] - 2019-08-19
### Added
- `_posture.csv` and `_angles.csv` are exported when `Export CSV` button is pressed.
  - `_posture.csv` table contains following five feature points (from head to tail direction)
    '  Columns: x1(m), x2(m), x3(m), x4(m), x5(m), y1(m), y2(m), y3(m), y4(m), y5(m)',
    '    x1(m): x coordinate of the head tip @ grid frame of reference',
    '    x2(m): x coordinate of the head-mid section @ grid frame of reference',
    '    x3(m): x coordinate of the mid section @ grid frame of reference',
    '    x4(m): x coordinate of the mid-tail section @ grid frame of reference',
    '    x5(m): x coordinate of the tail tip @ grid frame of reference',
    '    y1(m): y coordinate of the head tip @ grid frame of reference',        
  - `_angles.csv` table contains following angles
    '  Columns: a_hm(deg), a_tm(deg), a_bb(deg), a_tb(deg)', 
    '    a_hm(deg): head-mid section orientation (head half of the fish)', 
    '    a_tm(deg): tail-mid section orientation (tail half of the fish)', 
    '    a_bb(deg): body bend angle', 
    '    a_tb(deg): tail bend angle'

## [0.3.9] - 2018-11-22
### Added
- `Encounter` button is added which plots the trajectory with color-coded time since the last encounter

## [0.3.8] - 2018-10-15
### Changed
- change_log.txt renamed to changelog.md

### Fixed
- .cvs renamed to .csv (comma seperated value)
- trialset csv export tool warns if shape positions field does not exist.


## [0.3.7] - 2018-08-29
### Added
- "Learning curve" button exports .csv file containing path-length and durations
  - Dimensions: nSessions x nTrials
- "Bar plots" button exports .csv file containing stats by shapes for the probe trials
  - Rows: "Triangle Lg", "Triangle Sm", "Square Lg", "Square Sm", "Circle Lg", "Circle Sm", "Food", "Wall", "Not Wall", 
  - Columns: "Sampling density (counts/m)", "Sampling rate (Hz)", "Speed (m/s)", "Escan density (counts/m)", "Visit count", "Visit duration (s)", "Duration per visit", "Freq. backward swimming", 


## [0.3.6] - 2018-08-14
### Fixed
- "Import_Track" button works even when no _Track.mat files are found.
- "List files" button works even when no _Track.mat files are found.


## [0.3.5] - 2018-08-13
### Added
- "Bar plots"
  - Includes bar plots for the pooled probe trials
  - "csAnimals" is used to correctly select the animals to plot. (applies to List files, Learning curve)
- "Input Coordinates" button
  - Pressing "S" plots sampling density bar plots for the individual trial.
- Export CSV" button
  - Added velocity (signed) and distance/escan columns at the end
- "default.cfg"
  - Added distance to the circular wall from the center ("diameter_cm_wall")
    and exclusion distance ("dist_cm_wall") parameters


## [0.3.4] - 2018-08-12
### Added
- "Export to CSV (all)" and "Expport to CSV"
  - A_E (error angle: food vec-heading vec) now ranges from -90..90 deg (previously 0..90 deg)
- Bugfix: "Input coordinates"
  - [D]ensity keyboard shortcut added: plots the sampling density plot. 
  - [G]oto trial command now accepts a full file path as well as the trial ID.
  - Grid markers changed from '+' to '.' to permit right-clicking near the grid markers.
  - Closing the figure prompts the file saving dialog. Previously user had to press "OK" on the message box to save.

### Fixed
- trialset overview
  - Mouse click on the pixel now correctly copies the probe trial file location.
- "Check sync" button
  - "Shift" left/right jumps by 10 video frames and "Control" left/rigth jumps by 100 frames.


## [0.3.3] - 2018-08-10
### Added
- "Check Sync" button
- "Page Up/Down" button steps through the event markers.
- "Check Sync (all)" button 
  - Clicking the iamge cell displays the location of file and copies to the clipboard.
- "Sync timing" button
  - It applies the new time sync methods by using all the LED pulses (intead of using the first and last only).
- "Initialize" button
  - Faster video loading using the cached uncompressed video (_mov000x000.bin)
- "Import _Track" button
  - Copies _Track.mat files from one folder to each individual folder where the video file is located.
- "Export to CSV (all)" button
  - Data is resampled at 100 Hz after applying ADC-Camera time synchronization
    (previously sampled at each camera frame capture time)


## [0.3.2] - 2018-08-07
### Fixed
- "Fix Sync" button 
  - Manual time sync validation tool is added.
  - Outputs uncompressed movie preview (_mov000x000.bin) for faster subsequent loading.
- Trialset "Fix sync" button
   - Automatically corrects all sync errors in the trialset.

### Fixed
- "List files" button


## [0.3.1] - 2018-08-06
### Added
- "CSV export" feature now includes "_relations.csv" file that exports following info
  - Frame capture time (s)
  - distance to food (m)
  - Heading angle error (0..90 deg)
  - Proximity to each shape (boolean, controlled by "dist_cm_shapes" in defatul.cfg).
- New features in "Input Coordinates" tool
  - Counterclockwise rotation is added to the mouse right-click menu.
  - Enhanced image contrast.
  - Plot trajectory key is added (press "T" to toggle on/off).


## [0.3.0] - 2018-08-02
### Fixed
- Fixed csv export error. Orientation is correctly saved (previously divided by 20)

### Added
- Added a shortcut for "Input coordinates" figure: "[E]xport CSV" to export individual trial.


## [0.2.9] - 2018-08-01
- Trialset "Export to CSV" exports trajectories and shape locations in the grid frame of reference.
  - Rotation and scale factor are corrected and saved in metric and degree units.
  - (0,0) position corresponds to the center of the aquarium where (x,y) axis intersects.
  - Angle is measured from the x-axis in counterclockwise orientation.
  - Trajectory files are saved to "_Track.csv" and shape locations are saved to "_shapes.csv" file.


# [0.2.8] - 2018-07-31
### Fixed
- "Input Coordinates" oversize issue is corrected.
 

## [0.2.7] - 2018-07-31
### Added
- GUI "Density" button shows time spent per pixel
- "Learning curve" plot
  - "Input Coordinates" button allows a user to input the shape locations and orientatoins.
  - "Check FPS" button shows which trials have sync problem (fixable by loading the trial and hitting "Fix Sync" button).


## [0.2.6] - 2018-07-27
### Added
- GUI "Export" button: exporting trajectory to .csv files
  - Columns contain time(s), x-pos(m), y-pos(m), head-orientation (deg)
- "Export to CSV" in Trialset which exports the trajectory of each files to .csv
- "Probe Trials" in Trialset which collects shape coordinates and saves to "_Trialset.mat" file.


## [0.2.5] - 2018-07-26
### Fixed
- Fixed quantile_ error (now requires Statistics and Machine learning toolbox)
- Github update automatically overwrites changed files


## [0.2.4] - 2018-07-26
Compatible with Matlab R2014a version
Added FPS plot in Learning curve


## [0.2.3] - 2018-07-25
### Added
- GUI "Density" button made faster and improved
- Trialset "Bar plots" button added

## 0.2.2 - 2018-07-24
### Fixed
- Fixed missing settings_vistrack.m flie after software update.


## 0.2.1 - 2018-07-24
### Added
- .gitignore file added, ignoring *.wmv files


## 0.2.0 - 2018-07-24
### Fixed
- Debugged "visrack download-sample" command
  - Retrying up to 5 times


## 0.1.9 - 2018-07-24
### Fixed
- Debugged "visrack download-sample" command


## 0.1.8 - 2018-07-24
### Added
- GUI "Replay" window:
  - Added "[G]oto frame" keyboard shortcut.
  - Added fast frame jump by pressing "Shift" + Left/Right keys.
- Debugged "vistrack update" command.
- Uploaded sample files (R12A2*).
- Removed dependencies for "Statistics and Machine Learning" toolbox.
  - Only "Image processing toolbox" is required. 
  - "Parallel processing toolbox" is optional for faster processing.
- GUI "Initialize" button caches previous answers and a movie clip.
- New command added: "vistrack download-sample"
  - Download a sample video file that is too large to host from Github.com


## 0.1.7 - 2018-07-23
### Added
- Added "Summary" button in GUI
- Added "Export" button in GUI
- Robust video loading, number of frames calculated 3 times and max is taken.
- Updated "List Files" button for Trialset
- Added "Learning Curve" button for Trialset
- Settings box now displays new lines
- Faster video loading (~2x, requires a "Distrib_Computing_Toolbox")


## 0.1.6 - 2018-07-12
### Added
- Added .trialset file loading
- Added learing curve plot


## 0.1.5 - 2018-07-12
### Fixed
- Timestamp error is fixed when no ADC file is present.
- Initialize button error is fixed when the very first or last frame of the movie is selected.


## 0.1.4 - 2018-07-12
### Fixed
- Fixed update error
- Fixed tracking error


## 0.1.3 - 2018-07-12
### Added
- xyLED (LED coordinate) is detected automatically
- ADC recording doesn't need to be present
- Initialize button uses only the video, not ADC recording


## 0.1.2 - 2018-07-11
### Added
- vistrack.m script added
- Enabled save button when loading the tracking result


## 0.1.1 - 2018-07-10
### Added
- Initial commit