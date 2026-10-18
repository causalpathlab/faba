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
            0xc076570e130ab2a3,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xda2d28758bb10a9e,
            0xc0764e60fe4aced2,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xe93b7c4b3d4d76ba,
            0xc07556034cb5584d,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0x3d5e6ff958e9e321,
            0xc075a31eae77f04c,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xcaf522c8869429ff,
            0xc07625e9a9ad077e,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xacf3e0342f5c607a,
            0xc076b14716022419,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
    ],
    [
        [
            0x70e73ac83c7fbe47,
            0xc07523d3808ac864,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xcdcf9f1f7e376be0,
            0xc074c0ad5dec1b8e,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xcc02312a045e52b6,
            0xc0752dbb27c856d9,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0xbb9465a919de3be9,
            0xc0738ab91c730f21,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0x2ef2d51c0d21fd49,
            0xc073dde3060354cb,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
        [
            0x4c64cae062abcd14,
            0xc074b3893d872d05,
            0xd94d12186c0f2fb7,
            0x042a3110292eaa1d,
        ],
    ],
    [
        [
            0x1a6f86cb0258b91a,
            0xc08087635b88bbbd,
            0x21a0a17c74ee66e7,
            0xb1f2a56dffc98a7b,
        ],
        [
            0xead888dc14cf10ca,
            0xc0817b33788d58d0,
            0xc03fd980580c0016,
            0x729c211893bba020,
        ],
        [
            0x7d0e83a73ab2f607,
            0xc0819258f2264d42,
            0xe49eb533408eb935,
            0xe56c4e9227e166e2,
        ],
        [
            0x1c9051aad28c4fcf,
            0xc080ec33e98abcfa,
            0xb2909b8d4c1b4204,
            0x1de5a5a5311bada1,
        ],
        [
            0x4831d1dc59134a8d,
            0xc08228a97dc7e091,
            0x7455b4d6346fc030,
            0x83586eba2601d3e4,
        ],
        [
            0x500b2271ca1cb95c,
            0xc081d71b253aa22d,
            0x661e9550c5e0306c,
            0x001696de377771aa,
        ],
    ],
    [
        [
            0x1058faa8f21efb42,
            0xc07825bfc467e909,
            0x55227e6d270cad35,
            0x042a3110292eaa1d,
        ],
        [
            0xec864c6d677d6a69,
            0xc07897eab553e3c9,
            0x8960d0cd9d644356,
            0x042a3110292eaa1d,
        ],
        [
            0x6bb6a91f67330ac2,
            0xc078607498bdbe84,
            0xfdfeae043ad91c37,
            0x042a3110292eaa1d,
        ],
        [
            0xc45e8ddcc33648d4,
            0xc0781532c476f7ee,
            0x97266fbdfedcd013,
            0x042a3110292eaa1d,
        ],
        [
            0x478fef18f2287f47,
            0xc077af4bfe8886e2,
            0x7df181242cd6a259,
            0x042a3110292eaa1d,
        ],
        [
            0xa8effa32e9302aa1,
            0xc077ed722aa78a42,
            0x151ab972ed08ca72,
            0x042a3110292eaa1d,
        ],
    ],
];

/// Each fixture's chain is, bit for bit, the one recorded: its labels, log
/// joint, shifts and levels after every sweep. The bits are the platform
/// libm's (`exp`, `ln`, `ln_1p` are not correctly rounded everywhere), so
/// they hold where they were recorded, on Apple silicon, and elsewhere the
/// test is ignored.
#[test]
#[cfg_attr(
    not(all(target_os = "macos", target_arch = "aarch64")),
    ignore = "fingerprints recorded with the Apple silicon libm"
)]
fn fixed_chains_are_bit_for_bit_those_recorded() {
    for ((name, got), want) in chains().into_iter().zip(&RECORDED) {
        for (sweep, (g, w)) in got.iter().zip(want).enumerate() {
            assert_eq!(g, w, "{name} chain, sweep {sweep}: all {got:x?}");
        }
    }
}
