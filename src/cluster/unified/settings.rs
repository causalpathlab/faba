//! The blockmodel's settings, and what a run's record says of them.

/// The global weights over the cell roles.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Weights {
    /// `beta ~ Dirichlet(gamma / A)` with this `gamma`, drawn with the
    /// labels (a weak-limit hierarchical Dirichlet process).
    Hdp(f64),
    /// `beta = 1 / A`, fixed.
    Symmetric,
}

impl Weights {
    /// The HDP's `gamma`; none for symmetric weights.
    pub fn gamma(self) -> Option<f64> {
        match self {
            Weights::Hdp(g) => Some(g),
            Weights::Symmetric => None,
        }
    }
}

/// The model's sizes, concentrations and switches.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Settings {
    /// Cell roles `A`.
    pub cell_roles: usize,
    /// Gene roles `B`.
    pub gene_roles: usize,
    /// Dirichlet concentration of each cell's role weights around `beta`.
    pub alpha: f64,
    /// Dirichlet concentration of each gene's role weights.
    pub gene_alpha: f64,
    /// The cell roles' global weights: drawn (HDP) or fixed.
    pub weights: Weights,
    /// Prior standard deviation of each block's shift; `None`, no shift.
    pub shift_sd: Option<f64>,
    /// Profile strength `S` (infinite: the profile is the shifted baseline).
    #[serde(with = "strength")]
    pub strength: f64,
    /// Prior standard deviation of each cell's own level, which the model
    /// has at infinite strength only.
    pub cell_sd: f64,
}

impl Settings {
    /// Refuses settings that make the fit meaningless: roles out of range
    /// (labels are `u16`), or gene roles, or infinite strength, without the
    /// shift, which leave the roles no way to the reads.
    pub fn check(&self) -> anyhow::Result<()> {
        let roles = 1..u16::MAX as usize;
        anyhow::ensure!(
            roles.contains(&self.cell_roles) && roles.contains(&self.gene_roles),
            "cell roles and gene roles (--cell-roles, --gene-roles) must be between 1 and {}",
            u16::MAX - 1
        );
        anyhow::ensure!(
            self.gene_roles == 1 || self.shift_sd.is_some(),
            "{} gene roles need the block shift (--shift on): a site belongs to one \
             gene, so without it a gene's role never reaches its reads",
            self.gene_roles
        );
        anyhow::ensure!(
            self.strength.is_finite() || self.shift_sd.is_some(),
            "infinite profile strength needs the block shift (--shift on): without \
             it every block's reads are Binomial at the site's pooled fraction, so \
             no role reaches them"
        );
        Ok(())
    }
}

/// What a run's record says of the model it fitted (under `model` in its
/// results): its settings, and their preset when they are one (by its
/// `--model` name).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Recorded {
    pub preset: Option<String>,
    pub settings: Settings,
}

impl Recorded {
    /// In a few words: its preset's name, or that the settings are custom.
    pub fn label(&self) -> String {
        match &self.preset {
            Some(p) => format!("{p} preset"),
            None => "custom settings".to_string(),
        }
    }
}

/// A profile strength in JSON: a number, or `"inf"` (JSON has no infinity).
pub mod strength {
    use serde::{Deserialize, Deserializer, Serializer};

    /// `x` as a number, or `"inf"`.
    pub fn serialize<S: Serializer>(x: &f64, s: S) -> Result<S::Ok, S::Error> {
        if x.is_finite() {
            s.serialize_f64(*x)
        } else {
            s.serialize_str(&x.to_string())
        }
    }

    /// A number, or one written as text (`"inf"`).
    pub fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<f64, D::Error> {
        /// Either form.
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum Written {
            Number(f64),
            Text(String),
        }
        match Written::deserialize(d)? {
            Written::Number(x) => Ok(x),
            Written::Text(t) => t.parse().map_err(serde::de::Error::custom),
        }
    }
}
