coalesce_drop <- function(df){
  # Automatically coalesce all variables that end with .x and .y together and 
  # drops the .x and .y variants
  joined_vars <- df %>%
    select(ends_with(".x")) %>%
    names() %>%
    str_extract(".*(?=.x)")
  
  names(joined_vars) <- joined_vars
  
  df_new <- joined_vars %>%
    
    map_dfc(~ coalesce(
      df[[paste0(.x, ".x")]],
      df[[paste0(.x, ".y")]]
    ))
  
  df_final <- bind_cols(df, df_new)
  
  df_final <- df_final %>%
    select(-ends_with(c(".x", ".y"))) %>%
    distinct()
  return(df_final)
}
  