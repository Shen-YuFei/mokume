use std::{fmt, str::FromStr};

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum QuantMethod {
    DirectLfq,
    Ibaq,
    MaxLfq,
    TopN,
    Sum,
    Median,
    Ratio,
    Abd,
    Intensity,
    SpectralCount,
}

impl QuantMethod {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::DirectLfq => "directlfq",
            Self::Ibaq => "ibaq",
            Self::MaxLfq => "maxlfq",
            Self::TopN => "topn",
            Self::Sum => "sum",
            Self::Median => "median",
            Self::Ratio => "ratio",
            Self::Abd => "abd",
            Self::Intensity => "intensity",
            Self::SpectralCount => "spectral_count",
        }
    }
}

impl fmt::Display for QuantMethod {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

impl FromStr for QuantMethod {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "directlfq" => Ok(Self::DirectLfq),
            "ibaq" => Ok(Self::Ibaq),
            "maxlfq" => Ok(Self::MaxLfq),
            "topn" | "top3" => Ok(Self::TopN),
            "sum" | "all" | "allpeptides" => Ok(Self::Sum),
            "median" => Ok(Self::Median),
            "ratio" => Ok(Self::Ratio),
            "abd" | "abundance" | "tmtabundance" => Ok(Self::Abd),
            "intensity" | "reporter" | "tmtreporterintensity" => Ok(Self::Intensity),
            "spectral_count" | "spectralcount" | "count" => Ok(Self::SpectralCount),
            other => Err(format!("unknown quantification method: {other}")),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::QuantMethod;

    #[test]
    fn parses_aliases() {
        assert_eq!(
            QuantMethod::from_str("directlfq").ok(),
            Some(QuantMethod::DirectLfq)
        );
        assert_eq!(QuantMethod::from_str("all").ok(), Some(QuantMethod::Sum));
        assert_eq!(
            QuantMethod::from_str("spectralcount").ok(),
            Some(QuantMethod::SpectralCount)
        );
    }
}
