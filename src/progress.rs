use indicatif::{ProgressBar, ProgressStyle};

/// Create a progress bar of a specified length with desired styling.
pub fn progress_bar(len: u64) -> ProgressBar {
    let progress_bar = ProgressBar::new(len);
    progress_bar.set_style(ProgressStyle::default_bar().template(
        "[{elapsed_precise}] {bar:40.cyan/blue} {percent}% [{pos:>2}/{len:2}] [Remaining: {eta}]",
    ).expect("Invalid progress style."));

    progress_bar
}

/// Create a progress bar of a specified length and styling, with a terminal message.
pub fn progress_bar_msg(len: u64) -> ProgressBar {
    let progress_bar = ProgressBar::new(len);
    progress_bar.set_style(ProgressStyle::default_bar().template(
        "[{elapsed_precise}] {bar:20.cyan/blue} {percent}% [{pos:>2}/{len:2}] [Remaining: {eta}] [{msg}]",
    ).expect("Invalid progress style."));

    progress_bar
}
