fn main() {
    #[cfg(feature = "asm")]
    if std::env::consts::ARCH != "x86_64" {
        eprintln!("Currently feature `asm` can only be enabled on x86_64 arch.");
        std::process::exit(1);
    }
}
