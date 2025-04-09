use std::fs::File;
use std::path::PathBuf;

use anyhow::Context;
use log::{debug, LevelFilter};
use log4rs::append::console::{ConsoleAppender, Target};
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Root};
use log4rs::encode::pattern::PatternEncoder;
use log4rs::filter::threshold::ThresholdFilter;
use log4rs::{Config, Handle};
use tracing_appender::non_blocking::WorkerGuard;
use tracing_subscriber::fmt::layer;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{layer::SubscriberExt, Layer};

pub fn init_logging_smart(
    log_fp: Option<&PathBuf>,
    quiet_stdout: bool,
) -> Handle {
    let level = LevelFilter::Info;

    let file_endcoder = Box::new(PatternEncoder::new(
        "[{f}::{L}][{d(%Y-%m-%d %H:%M:%S)}][{l}] {m}{n}",
    ));
    let console_encoder = Box::new(PatternEncoder::new("{h(>)} {m}{n}"));
    let stderr = ConsoleAppender::builder()
        .encoder(console_encoder)
        .target(Target::Stderr)
        .build();

    let config = if let Some(fp) = log_fp {
        let logfile = FileAppender::builder()
            .encoder(file_endcoder)
            .build(fp)
            .expect("failed to make file logger");
        let mut config = Config::builder();
        let logfile_appender =
            Appender::builder().build("logfile", Box::new(logfile));
        config = config.appender(logfile_appender);
        let mut root_logger = Root::builder().appender("logfile");
        if !quiet_stdout {
            config = config.appender(
                Appender::builder()
                    .filter(Box::new(ThresholdFilter::new(level)))
                    .build("stderr", Box::new(stderr)),
            );
            root_logger = root_logger.appender("stderr");
        }

        config.build(root_logger.build(LevelFilter::Trace)).unwrap()
    } else if !quiet_stdout {
        Config::builder()
            .appender(
                Appender::builder()
                    .filter(Box::new(ThresholdFilter::new(level)))
                    .build("stderr", Box::new(stderr)),
            )
            .build(Root::builder().appender("stderr").build(LevelFilter::Trace))
            .unwrap()
    } else {
        Config::builder()
            .build(Root::builder().build(LevelFilter::Trace))
            .unwrap()
    };

    let handle = log4rs::init_config(config).expect("failed to init logging");
    let command_line = std::env::args().collect::<Vec<String>>().join(" ");
    debug!("command line: {command_line}");
    handle
}

pub fn init_logging(log_fp: Option<&PathBuf>) -> Handle {
    init_logging_smart(log_fp, false)
}

pub fn init_tracing(
    log_fp: Option<&PathBuf>,
) -> anyhow::Result<Option<WorkerGuard>> {
    let mut layers = Vec::new();
    let stderr_layer = layer()
        .compact()
        .without_time()
        .with_target(false)
        .with_filter(tracing::level_filters::LevelFilter::INFO)
        .boxed();
    layers.push(stderr_layer);

    let guard = if let Some(fp) = log_fp {
        let log_fh =
            File::create(fp).context("failed to create output file")?;
        let (logfile_appender, guard) = tracing_appender::non_blocking(log_fh);

        let file_layer = layer()
            .json()
            .with_target(true)
            .with_line_number(true)
            .with_file(true)
            .with_writer(logfile_appender)
            .with_filter(tracing::level_filters::LevelFilter::DEBUG)
            .boxed();
        layers.push(file_layer);
        Some(guard)
    } else {
        None
    };
    tracing_subscriber::registry().with(layers).init();
    let command_line = std::env::args().collect::<Vec<String>>().join(" ");
    debug!("command line: {command_line}");
    Ok(guard)
}
