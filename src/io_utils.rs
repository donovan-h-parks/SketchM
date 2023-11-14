use std::ffi::{OsStr, OsString};
use std::path::{Path, PathBuf};

/// Append extension to path if it doesn't already end with the extension.
pub fn append_extension_to_path(path: &Path, extension: &str) -> PathBuf {
    let mut path_ext = path.to_path_buf();

    if path.extension() != Some(OsStr::new(extension)) {
        let mut out_file_ext: OsString = path.into();
        out_file_ext.push(".");
        out_file_ext.push(extension);
        path_ext = out_file_ext.into();
    }

    path_ext
}
