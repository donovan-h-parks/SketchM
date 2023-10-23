use log::LevelFilter;
use log4rs::{
    append::console::{ConsoleAppender, Target},
    config::{Appender, Config, Root},
    encode::pattern::PatternEncoder,
    filter::threshold::ThresholdFilter,
};

/// Configure logger to write to stderr.
pub fn setup_logger() {
    let level = log::LevelFilter::Info;
    let pattern = "[{d(%Y-%m-%d %H:%M:%S)}] {h({l})}: {m}{n}";

    // log to stderr
    let stderr = ConsoleAppender::builder()
        .encoder(Box::new(PatternEncoder::new(pattern)))
        .target(Target::Stderr)
        .build();

    // configure logging
    let config = Config::builder()
        .appender(
            Appender::builder()
                .filter(Box::new(ThresholdFilter::new(level)))
                .build("stderr", Box::new(stderr)),
        )
        .build(Root::builder().appender("stderr").build(LevelFilter::Trace))
        .expect("Failed to configure logger.");

    log4rs::init_config(config).expect("Failed to initialize logger.");
}
