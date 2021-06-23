############################################
### TraceLab Analysis Configuration File ###
############################################

### DEMI-Specific Settings ###

# Whether or not interpolated channels should be dropped
drop_interpolated <- FALSE

# Whether or not EEG data has been Current Source Density (CSD) processed
using_csd <- TRUE

# The frequency bins (in Hz) to use for wavelet-decomposition of the EEG signal
wt_frequencies <- 1:40

# Whether to use the 4 seconds before each trial's onset as an EEG baseline
# instead of the windows 500 ms to 200 ms before each epoch
use_pretrial_baseline <- TRUE



### General Settings ###

# The resolution of the screen data was collected with.
screen_res <- c(1920, 1080)

# Whether to render plots of points/trials caught by tracing filters
# (useful for initial setup/adjustment, but otherwise slows down
# preprocessing a ton)
plot_filters <- TRUE



### Tracing Filter Settings ###

## 'No Shape' Filter

# On some trials, the tracing consists only of a few samples with no clear
# shape (either a short line or a small cluster near the origin). This filter
# tries to catch and drop these using the figure-to-tracing size ratio.
#
# - 'min_size_ratio': The minimum ratio of figure size to tracing size for the
#    tracing to be flagged as a "no shape" trial. For this filter, "size"
#    refers to whichever dimension is largest (height or width).

no_shape_params <- list(
  min_size_ratio = 5
)


## "Trial Done" Filter

# On some trials, participants attempt to end the trial but miss the origin,
# so they either loop back to it or lift their finger and re-attempt. This
# filter tries to detect when this happens and drop all points after the
# failed "end trial" attempt from further analysis.
#
# - 'origin_radius': An expanded radius around the origin (in px). If the
#    tracing passes within this radius and the proportion of tracing complete
#    is >= 'min_prop', the tracing will be considered done.
#
# - 'end_radius': A radius around the origin (in px) all remaining points must
#    be within before the tracing can be considered done.
#
# - 'pause_radius': The maximum distance (in px) that the previous point can be
#    from the origin for a sufficiently long pause (as defined by 'min_pause')
#    to be considered a tracing end.
#
# - 'min_prop': The minimum proportion of the tracing that needs to be complete
#    before the tracing can be considered done.
#
# - 'end_prop': The minimum proportion of the tracing that needs to be complete
#    before a sufficiently long pause can be considered a trial end.
#
# - 'min_pause': The minimum time difference (in ms) between two points in the
#    end proportion of the trace for the trial to be considered done (provided
#    that the previous sample is within 'pause_radius' of the origin).

done_filter_params <- list(
  origin_radius = 50,
  end_radius = 400,
  pause_radius = 200,
  min_prop = 0.62,
  end_prop = 0.82,
  min_pause = 0.100
)

# For this dataset, we do a second pass of this filter with different parameters

done_filter_params2 <- list(
  origin_radius = 50,
  end_radius = 400,
  pause_radius = 200,
  min_prop = 0.62,
  end_prop = 0.55,
  min_pause = 0.400
)



## Touchscreen Glitch Filter

# Occasionally during a trial, the touchscreen glitches out and jumps the cursor
# to a seemingly random location on the screen (can occur if the participant
# lifts their finger at all during the tracing). This filter tries to catch and
# remove these glitch points so those trials don't need to be excluded.
#
# This filter has four parts: a filter that detects glitches based on angle and
# distance, two filters that detect whether the first or last points of the
# tracing are glitch points, respectively, and a final filter that detects
# multiple glitches in a row when they share x or y positions.
#
# - 'min_angle_diff': The minimum absolute change in direction (in degrees) from
#    the jump towards a point vs the jump away from it for a point to be able to
#    be flagged as an angle glitch.
#
# - 'min_angle_diff_alt': The minimum absolute change in direction (in degrees)
#    from the jump towards a point vs the jump to the point *two after* (i.e.
#    ignoring a subsequent glitch point) for an angle_glitch to be considered
#    valid. For reducing false positives.
#
# - 'min_dist': The minimum size (in px) that the jumps towards and away from a
#    point must be for a point to be able to be flagged as an angle glitch. Also
#    used for start glitches (min. dist first point must be away from origin and
#    min. dist first line segment needs to be for first point to be glitch) and
#    end glitches (min. increase in dist from origin vs 2nd-last point for end
#    point to be glitch).

glitch_filter_params <- list(
  min_angle_diff = 88,
  min_angle_diff_alt = 70,
  min_dist = 173
)


## False Start Filter

# Sometimes participants have issues getting a trial to start, which can lead
# to some small weird clusters of samples separated from the actual tracing
# start by a small-to-substantial time delay (and can cause otherwise-good
# trials to fail the 'excessive gap' filter below if the time delay is large
# enough). This filter checks for time jumps within the first part of the
# tracing and drops all samples prior to the last eligible jump.
#
# - 'start_radius': The maximum distance from the origin that the first "real"
#    sample of the tracing can have for previous points to be flagged as "false
#    start".
#
# - 'max_prop': The maximum proportion of the tracing that can be complete for
#    prior samples to be flagged as "false start".
#
# - 'min_pause': The minimum time difference (in ms) between points for
#    previous points to be dropped.

false_start_params <- list(
  start_radius = 80,
  max_prop = 0.20,
  min_pause = 0.064
)


## Hand Noise Filter

# During a tracing, a participant's palm or forearm can sometimes briefly make
# contact with the screen, causing a gap and a small group of samples off of the
# normal tracing path. This filter attempts to flag hand noise when it occurs,
# potentially allowing the trial to be saved if the gap is small enough.
#
# This filter consists of two parts: an initial filter looking for at least one
# above-threshold large & sharp jump between points during each trial, and a
# subsequent filter to look for and flag hand noise on those trials only. The
# latter filter is based on flagging all points between the two largest jumps
# in a tracing, provided that those jumps meet certain thresholds. If there are
# multiple instances of hand noise in a trial, this filter will not work
# properly.
#
# - 'min_sharp_turnsum': The minimum sum of absolute change in direction
#    (in degrees) of the current and next points for a point to be considered
#    a "sharp jump", provided that the 'min_dist_a' threshold is also met
#    for the same point.
#
# - 'max_size': The maximum allowed height or width of a detected hand noise
#    region. Flagged regions larger than this size will be ignored.
#
# - 'max_timediff': The maximum time gap (in ms) that the two largest jumps in
#    a tracing can have. If either jump takes longer than this, the points
#    between them will not be flagged as hand noise.
#
# - 'min_dist_a': The minimum size (in px) for each of the two largest jumps.
#
# - 'min_dist_b': The minimum size (in px) for each of the two largest jumps
#    if the direction change of the jump is below the 'min_angle_diff'
#    threshold.
#
# - 'min_angle_diff': The minimum absolute change in direction (in degrees)
#    for each of the two largest jumps if the distance of the jump is above
#    the 'min_dist_a' threshold but below the 'min_dist_b' one.
#
# - 'min_end_dist': The minimum size (in px) of jump away from the trial origin
#    for a cluster of points at the end of a trial to be flagged as hand noise.

hand_noise_params <- list(
  min_sharp_turnsum = 90,
  max_size = 100,
  max_timediff = 0.100,
  min_dist_a = 120,
  min_dist_b = 300,
  min_angle_diff = 30,
  min_end_dist = 100
)


## Incomplete Trial Filter

# On some trials, participants accidentally end their tracing early by tracing
# too close to the origin point ("accidental end" trials). On others, the
# touchscreen loses track of their finger part-way, leaving an incomplete shape
# ("end gap" trials). This filter tries to catch and exclude both kinds of
# trials from further analysis.
#
# "Accidental end" trials are filtered out using either the 'min_size_ratio'
# parameter or a combination of the 'min_length_ratio' and 'min_shift_diff'
# parameters. "End gap" trials are filtered out using the 'min_end_gap'
# parameter alone.
#
# - 'min_end_gap': The minimum size (in px) of the gap between the first and
#    last points of a tracing for it to be flagged as an "end gap" trial.
#
# - 'min_size_ratio': The minimum ratio of figure size to tracing size for the
#    tracing to be flagged as an "accidental end" trial. For this filter,
#    "size" refers to whichever dimension is largest (height or width).
#
# - 'min_length_ratio': The minimum ratio of figure path length to tracing path
#    length for a trial to be flagged as an "accidental end" trial (if the
#    'min_shift_diff' threshold is also met).
#
# - 'min_shift_diff': The minimum difference in absolute lateral shift between
#    the tracing and figure for a trial to be flagged as an "accidental end"
#    trial (if the 'min_length_ratio' threshold is also met). "Lateral shift"
#    here refers to how much a figure is shifted left or right from being
#    centered on the origin point: for example, a tracing where all samples are
#    to the left of the origin would have a lateral shift of >= 1, whereas a
#    a tracing centered directly on the origin would have a shift of 0. If a
#    figure is fairly centered on origin but its tracing is not, this is a
#    good indicator of an accidental trial end.

incomplete_params <- list(
  min_end_gap = 300,
  min_size_ratio = 2.5,
  min_length_ratio = 1.45,
  min_shift_diff = 0.71
)


## Excessive Gap Filter

# Gaps in tracings get interpolated after filtering, so this filter tries to
# flag trials with gaps too large (by time or by distance) to be safely
# interpolated.
#
# This filter has four levels to it:
#  a) A time-only filter for catching extra-long pauses
#  b) A time + distance filter for catching *really* large gaps
#  c) A time + distance + angle filter for catching large gaps during turns
#  d) A time + distance + angle filter for catching gaps during sharp turns
#
# - 'min_pause': The longest allowable pause (in ms) between two samples, beyond
#    which a trial will be flagged for exclusion.
#
# - 'min_timegap': The minimum time difference (in ms) between two samples for
#    any distance between them to be considered an unrecoverable gap.
#
# - 'min_dist_b': The minimum distance (in px) between two points for the gap
#    to be considered unrecoverable, regardless of angle.
#
# - 'min_dist_c': The minimum distance (in px) between two points for the gap
#    to be considered unrecoverable, provided that the cumulative absolute
#    angle change between the pre-gap and post-gap segments is above the
#    threshold set by 'min_turnsum_c'.
#
# - 'min_dist_d': The minimum distance (in px) between two points for the gap
#    to be considered unrecoverable, provided that the cumulative absolute
#    angle change between the pre-gap and post-gap segments is above the
#    threshold set by 'min_turnsum_d'.
#
# - 'min_turnsum_c': The minimum sum (in degrees) of absolute change in
#    direction of the pre and post-gap samples for the gap to be considered
#    unrecoverable, provided that the gap is larger than the threshold set
#    by 'min_dist_c'.
#
# - 'min_turnsum_d': The minimum sum (in degrees) of absolute change in
#    direction of the pre and post-gap samples for the gap to be considered
#    unrecoverable, provided that the gap is larger than the threshold set
#    by 'min_dist_d'.

gap_filter_params <- list(
  min_pause = 0.170,
  min_timegap = 0.070,
  min_dist_b = 500,
  min_dist_c = 250,
  min_dist_d = 180,
  min_turnsum_c = 90,
  min_turnsum_d = 120
)
