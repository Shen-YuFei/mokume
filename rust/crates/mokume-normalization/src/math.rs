use mokume_core::{MokumeError, Result};

pub(crate) fn valid_scale(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

pub(crate) fn unsupported<T>(stage: &'static str) -> Result<T> {
    Err(MokumeError::NotImplemented { stage })
}
