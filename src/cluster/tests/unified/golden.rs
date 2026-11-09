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
            0xc076571b13389fcb,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xda2d28758bb10a9e,
            0xc0764e6d90f91bc9,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xe93b7c4b3d4d76ba,
            0xc07556105153e70b,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0x3d5e6ff958e9e321,
            0xc075a32a8aee9721,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xcaf522c8869429ff,
            0xc07625f5bcdd2264,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xacf3e0342f5c607a,
            0xc076b14eef9e754e,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
    ],
    [
        [
            0x70e73ac83c7fbe47,
            0xc07523de241256b5,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xcdcf9f1f7e376be0,
            0xc074c0bdfd941dd0,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xcc02312a045e52b6,
            0xc0752dcb400b26e9,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xbb9465a919de3be9,
            0xc0738ac6c7510a92,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0x2ef2d51c0d21fd49,
            0xc073dded09890276,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0x4c64cae062abcd14,
            0xc074b39b0caf2d00,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
    ],
    [
        [
            0x1a6f86cb0258b91a,
            0xc080879beb82b40b,
            0x5c8135245f140d15,
            0x4d180ff6bdf15c9b,
        ],
        [
            0x6639799233dc7083,
            0xc080f8d738c71d01,
            0x4aee1efffa260598,
            0x44399ac6c9e28fb1,
        ],
        [
            0x1e791fc033b579d9,
            0xc08136b67133bb7e,
            0x290b52a96b92e2f2,
            0x8e90f2ad807f1137,
        ],
        [
            0xbc2c9a491b0450e3,
            0xc0816439370c5209,
            0x2a56e293b5cb051a,
            0xa2036219faea38fb,
        ],
        [
            0xe0259a0f0daa1503,
            0xc081778e6bf67a21,
            0x137df543aecd3518,
            0x59d8c53e29492ff7,
        ],
        [
            0x7a17ba7a09fe5355,
            0xc081c2c298e85bea,
            0x3f9d7a4a1b0e29a1,
            0xb4b0998787c016a5,
        ],
    ],
    [
        [
            0x1058faa8f21efb42,
            0xc07825c442d43381,
            0x55227e6d270cad35,
            0x042a3110292eaa1d,
        ],
        [
            0xec864c6d677d6a69,
            0xc07897ef95f199f5,
            0x8960d0cd9d644356,
            0x042a3110292eaa1d,
        ],
        [
            0x6bb6a91f67330ac2,
            0xc07860764d372f22,
            0xfdfeae043ad91c37,
            0x042a3110292eaa1d,
        ],
        [
            0xc45e8ddcc33648d4,
            0xc0781534d3232f11,
            0x97266fbdfedcd013,
            0x042a3110292eaa1d,
        ],
        [
            0x478fef18f2287f47,
            0xc077af56c8322d8f,
            0x7df181242cd6a259,
            0x042a3110292eaa1d,
        ],
        [
            0xa8effa32e9302aa1,
            0xc077ed7bb9dce6fd,
            0x151ab972ed08ca72,
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
