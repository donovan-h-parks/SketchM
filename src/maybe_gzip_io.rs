use std::io::{BufWriter, Write};
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::Result;
use gzp::{
    deflate,
    par::compress::{ParCompress, ParCompressBuilder},
    ZWriter,
};
use flate2::read::GzDecoder;


const BUFFER_SIZE: usize = 128 * 1024;

pub struct GzipParams {
    pub level: usize,
    pub threads: usize,
}

pub enum MaybeGzipWriter<T: Write> {
    Gzip(ParCompress<deflate::Gzip>),
    Plain(BufWriter<T>),
}

impl<T: Write + Send + 'static> MaybeGzipWriter<T> {
    pub fn new(file: T, gzip_params: Option<GzipParams>) -> Self {
        match gzip_params {
            Some(params) => MaybeGzipWriter::Gzip(
                ParCompressBuilder::new()
                    .compression_level(flate2::Compression::new(params.level as u32))
                    .num_threads(params.threads)
                    .expect("Invalid thread count")
                    .from_writer(file),
            ),
            None => MaybeGzipWriter::Plain(BufWriter::with_capacity(BUFFER_SIZE, file)),
        }
    }

    pub fn finish(self) -> Result<()> {
        match self {
            MaybeGzipWriter::Gzip(mut gzip) => gzip.finish()?,
            MaybeGzipWriter::Plain(mut bw) => bw.flush()?,
        };

        Ok(())
    }
}

impl<T: Write> MaybeGzipWriter<T> {
    fn as_write(&mut self) -> &mut dyn Write {
        match self {
            MaybeGzipWriter::Gzip(gzip) => gzip,
            MaybeGzipWriter::Plain(plain) => plain,
        }
    }
}

impl<T: Write> Write for MaybeGzipWriter<T> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.as_write().write(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.as_write().flush()
    }
}

/// Read normal or compressed files based on the absence
/// or presence of a `gz` extension.
pub fn maybe_gzip_reader(filename: &str) -> Result<Box<dyn BufRead + Send>> {
    let path = Path::new(filename);
    let file = File::open(path)?;

    if path.extension() == Some(OsStr::new("gz")) {
        Ok(Box::new(BufReader::with_capacity(
            BUFFER_SIZE,
            GzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(BUFFER_SIZE, file)))
    }
}
