//! What the read views draw of a site's reads: the converted ones (the
//! default), the unconverted ones, both (converted in front of the total),
//! or the sites themselves. `c` cycles them in `faba pileup` and `faba
//! metagene`.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Show {
    #[default]
    Converted,
    Unconverted,
    Both,
    Sites,
}

impl Show {
    pub fn next(self) -> Self {
        match self {
            Show::Converted => Show::Unconverted,
            Show::Unconverted => Show::Both,
            Show::Both => Show::Sites,
            Show::Sites => Show::Converted,
        }
    }

    /// The y axis, given the channels' names (`on` converted, `off` not).
    pub fn label(self, on: &str, off: &str) -> String {
        match self {
            Show::Converted => format!("{on} reads"),
            Show::Unconverted => format!("{off} reads"),
            Show::Both => format!("{on} in front of {on} + {off} reads"),
            Show::Sites => "sites".into(),
        }
    }

    /// A short name, for the key that switches to it.
    pub fn short(self) -> &'static str {
        match self {
            Show::Converted => "converted",
            Show::Unconverted => "unconverted",
            Show::Both => "both",
            Show::Sites => "sites",
        }
    }
}
