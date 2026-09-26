use super::*;
use crate::figure::{Anchor, Canvas, INK};
use ratatui::backend::TestBackend;
use ratatui::Terminal;

fn svg(w: f64, h: f64) -> String {
    let mut c = Canvas::new(w, h);
    c.text(4.0, 12.0, "PANEL1", 9.0, Anchor::Start, INK);
    c.finish()
}

#[test]
fn raster_has_the_scaled_size() {
    let img = raster(&svg(100.0, 40.0), 2.0).unwrap();
    assert_eq!((img.width(), img.height()), (200, 80));
}

#[test]
fn plot_image_rebuilds_only_when_stale_or_resized() {
    let picker = Picker::halfblocks();
    let mut plot = PlotImage::default();
    let mut term = Terminal::new(TestBackend::new(40, 10)).unwrap();
    let mut calls = 0;
    let mut frame = |plot: &mut PlotImage, area: Rect, calls: &mut i32| {
        term.draw(|f| {
            let drawn = plot.render(f, area, &picker, |w, h| {
                *calls += 1;
                svg(w, h)
            });
            assert!(drawn);
        })
        .unwrap();
    };
    let area = Rect::new(0, 0, 20, 5);
    frame(&mut plot, area, &mut calls);
    frame(&mut plot, area, &mut calls);
    assert_eq!(calls, 1, "unchanged view: cached");
    plot.invalidate();
    frame(&mut plot, area, &mut calls);
    assert_eq!(calls, 2, "stale: rebuilt");
    frame(&mut plot, Rect::new(0, 0, 30, 6), &mut calls);
    assert_eq!(calls, 3, "resized: rebuilt");
}

#[test]
fn only_known_image_terminals_are_queried() {
    let env = |pairs: &'static [(&'static str, &'static str)]| {
        move |k: &str| {
            pairs
                .iter()
                .find(|(key, _)| *key == k)
                .map(|(_, v)| v.to_string())
        }
    };
    assert!(likely_images(env(&[("TERM_PROGRAM", "iTerm.app")])));
    assert!(likely_images(env(&[("TERM", "xterm-kitty")])));
    assert!(likely_images(env(&[("KITTY_WINDOW_ID", "3")])));
    assert!(likely_images(env(&[
        ("TERM", "tmux-256color"),
        ("LC_TERMINAL", "iTerm2")
    ])));
    assert!(!likely_images(env(&[
        ("TERM_PROGRAM", "Apple_Terminal"),
        ("TERM", "xterm-256color")
    ])));
    assert!(!likely_images(env(&[("TERM", "screen-256color")])));
    assert!(!likely_images(env(&[("KITTY_WINDOW_ID", "")])));
}
