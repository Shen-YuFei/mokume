use crate::{MokumeError, Result};

pub fn parse_memory_to_bytes(value: &str) -> Result<u64> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        return Err(MokumeError::InvalidMemory {
            value: value.to_string(),
        });
    }

    let split_at = trimmed
        .find(|ch: char| !(ch.is_ascii_digit() || ch == '.'))
        .unwrap_or(trimmed.len());
    let number = &trimmed[..split_at];
    let unit = trimmed[split_at..].trim().to_ascii_uppercase();

    let quantity = number
        .parse::<f64>()
        .map_err(|_| MokumeError::InvalidMemory {
            value: value.to_string(),
        })?;
    if !quantity.is_finite() || quantity <= 0.0 {
        return Err(MokumeError::InvalidMemory {
            value: value.to_string(),
        });
    }

    let multiplier = match unit.as_str() {
        "" | "B" => 1_u64,
        "K" | "KB" => 1024_u64,
        "M" | "MB" => 1024_u64.pow(2),
        "G" | "GB" => 1024_u64.pow(3),
        "T" | "TB" => 1024_u64.pow(4),
        _ => {
            return Err(MokumeError::InvalidMemory {
                value: value.to_string(),
            });
        }
    };

    let bytes = quantity * multiplier as f64;
    if !bytes.is_finite() || bytes > u64::MAX as f64 || bytes.round() < 1.0 {
        return Err(MokumeError::InvalidMemory {
            value: value.to_string(),
        });
    }
    Ok(bytes.round() as u64)
}

pub fn parse_memory_to_gib(value: &str) -> Result<f64> {
    Ok(parse_memory_to_bytes(value)? as f64 / 1024_f64.powi(3))
}

#[cfg(test)]
mod tests {
    use super::{parse_memory_to_bytes, parse_memory_to_gib};

    #[test]
    fn parses_common_memory_values() {
        assert_eq!(parse_memory_to_bytes("1GB").ok(), Some(1024_u64.pow(3)));
        assert_eq!(
            parse_memory_to_bytes("512MB").ok(),
            Some(512 * 1024_u64.pow(2))
        );
        assert_eq!(parse_memory_to_gib("80GB").ok(), Some(80.0));
    }

    #[test]
    fn rejects_invalid_memory_values() {
        assert!(parse_memory_to_bytes("4XB").is_err());
        assert!(parse_memory_to_bytes("0GB").is_err());
        assert!(parse_memory_to_bytes("0.1B").is_err());
    }
}
