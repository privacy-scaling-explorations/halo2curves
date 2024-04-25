mod field;
mod utils;

#[proc_macro]
pub fn impl_field(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    field::impl_field(input)
}
