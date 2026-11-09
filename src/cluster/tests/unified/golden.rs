//! Fixed chains, fingerprinted sweep by sweep: a change to the sampler's
//! arithmetic that should leave its every bit alone is held to that, on the
//! platform the fingerprints were recorded on.

use super::*;

/// The fixtures' chains, each `(name, model's fingerprints)`: the link
/// preset over whole reads and over fractionally weighted ones, the blocks
/// preset, and finite strength with gene roles and the shift over weighted
/// reads (the shift's draws through the pseudo-counts).
fn chains() -> Vec<(&'static str, Vec<[u64; 4]>)> {
    let (reads, gene_of, _) = planted_links(30, 6, 171);
    let whole = Links::new(&Weighted::new(&reads, None), &gene_of);
    let mut rng = SmallRng::seed_from_u64(172);
    let mut weight =
        |n: usize| -> Vec<f64> { (0..n).map(|_| rng.random_range(0.1..1.9)).collect() };
    let (on, off) = (weight(reads.on.len()), weight(reads.off.len()));
    let weighted = Links::new(&Weighted::new(&reads, Some((&on, &off))), &gene_of);
    let (block_reads, block_genes, _, _, _) = planted_blocks(30, 6, 1.0, 173);
    let blocks = Links::new(&Weighted::new(&block_reads, None), &block_genes);
    let run = |links: &Links, genes: &[u32], s: Settings| {
        Model::new(links, genes, genes.len(), s)
            .expect("an allowed setting")
            .fingerprints(SWEEPS, 174)
    };
    vec![
        ("link", run(&whole, &gene_of, link_preset(3, 0.5))),
        (
            "weighted link",
            run(&weighted, &gene_of, link_preset(3, 0.5)),
        ),
        ("blocks", run(&blocks, &block_genes, blocks_preset())),
        (
            "shifted weighted",
            run(
                &weighted,
                &gene_of,
                Settings {
                    cell_roles: 3,
                    ..settings(Weights::Hdp(1.0), true, 2.0, 2)
                },
            ),
        ),
    ]
}

/// Sweeps fingerprinted per chain.
const SWEEPS: usize = 6;

/// The chains' fingerprints as recorded, in [`chains`]' order.
const RECORDED: [[[u64; 4]; SWEEPS]; 4] = [
    [
        [
            0x4780f97488948d9d,
            0xc07661e87db81fca,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xa4052cf0c82906f1,
            0xc07515f95089371f,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0x8e65543ad1a355f3,
            0xc076bbcd6a4970e8,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0x1a71c3ca287d6108,
            0xc077ce09a3216858,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0x8b8b5ce023a0ecf5,
            0xc0771abcd75ee3c7,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xa5a3f713c1a9b5d6,
            0xc0771971f519f92e,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
    ],
    [
        [
            0x70e73ac83c7fbe47,
            0xc0752121445a56b5,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0x421e3e327c7a4b7a,
            0xc074424d07d861cb,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xebd454ea25ac3b9e,
            0xc074615b56221ce1,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0x2bea07ea83b9a206,
            0xc074cbd1be889813,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xd6bae3109533f20c,
            0xc073ecd7d708101d,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0x24d66aada1950a6d,
            0xc0758b355dd85259,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
    ],
    [
        [
            0x1a6f86cb0258b91a,
            0xc0807f2ec7b1ce80,
            0x6f0d7f660950ae4f,
            0x523137655bf59827,
        ],
        [
            0x76e336107700b900,
            0xc081e289d38a285c,
            0xa02d7557fc82beca,
            0x416c5d565ef8c247,
        ],
        [
            0x5e9286c00caba282,
            0xc08171855e7841ad,
            0xc3def8386795391a,
            0xd5b7d7b200d50a32,
        ],
        [
            0x7bd28dee4f02ab77,
            0xc081b191ec6d9c7a,
            0xb2d72a2bd0505371,
            0x77f6d9b5205e1ea8,
        ],
        [
            0xb48c30d81c83204d,
            0xc08100903979c9e1,
            0x5712d8ef2c312f67,
            0x1b7613828a696829,
        ],
        [
            0x538a343015274e18,
            0xc080ed7a9b02ce05,
            0x76e29557d74d5220,
            0xc0a7de91b9922a79,
        ],
    ],
    [
        [
            0x1058faa8f21efb42,
            0xc078510042c011e8,
            0xfa16e6c649ceb895,
            0x042a3110292eaa1d,
        ],
        [
            0x61bcc91bafc0b367,
            0xc078610dfb34e325,
            0x7e6fce503d03ced0,
            0x042a3110292eaa1d,
        ],
        [
            0xf31e5ac1cc05d53f,
            0xc077fcffd2ef9199,
            0x43de8eb0aebfcbd2,
            0x042a3110292eaa1d,
        ],
        [
            0x7e45547b46e82daa,
            0xc0777a85172c0218,
            0x2fb529e9b5c7ef0b,
            0x042a3110292eaa1d,
        ],
        [
            0x212174530adfad92,
            0xc075e2071003b29f,
            0x234fc6a8a906fdfc,
            0x042a3110292eaa1d,
        ],
        [
            0x89e5ff947217d671,
            0xc07789d2e9b15a1a,
            0x5256a33063e0d6cd,
            0x042a3110292eaa1d,
        ],
    ],
];

/// Each fixture's chain is, bit for bit, the one recorded: its labels, log
/// joint, shifts and levels after every sweep. The bits are the platform
/// libm's (`exp`, `ln`, `ln_1p` are not correctly rounded everywhere), so
/// they hold where they were recorded, on Linux x86_64, and elsewhere the
/// test is ignored.
#[test]
#[cfg_attr(
    not(all(target_os = "linux", target_arch = "x86_64")),
    ignore = "fingerprints recorded with the Linux x86_64 libm"
)]
fn fixed_chains_are_bit_for_bit_those_recorded() {
    for ((name, got), want) in chains().into_iter().zip(&RECORDED) {
        for (sweep, (g, w)) in got.iter().zip(want).enumerate() {
            assert_eq!(g, w, "{name} chain, sweep {sweep}: all {got:x?}");
        }
    }
}
