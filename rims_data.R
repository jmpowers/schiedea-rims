library(tidyverse)
library(googlesheets4)
gsheet <- gs4_get("1pYbAnEDw2KfM34l85wlJV6pfAr1DroPj_7GjfApnCq8")
sheet_names(gsheet)

hk.pops <- set_names(c("WK","WK","WK","879WKG","892WKG","904WPG","3587WP"),
                     c("794","866","899","879","892","904","3587")) #merge Waianae Kai populations (WK)

add_combos <- function(x) {
  mutate(x, 
         sxc = paste(species, crosstype, sep="."),
         sxcxm = paste(species, crosstype, mompop, momid, sep="."),
         mompid = paste(mompop, momid, sep="."),
         across(any_of("dadid"), dadpid= ~paste(dadpop, ., sep=".")),
         smompop = paste(species, mompop, sep=""),
         across(where(is.character), as.factor),
         crosstype = relevel(crosstype, "hybrid", 3),
         across(any_of(c("mompop", "dadpop")), recode_factor, !!!hk.pops)) 
}

# hybrids -----------------------------------------------------------------

hybrids <- read_sheet(gsheet, "hybrids", col_types="c")

# seeds -------------------------------------------------------------------

seeds <- read_sheet(gsheet, "seeds", col_types="c") %>% 
  filter(dadpop != "closed") %>% 
  mutate(seednum = as.integer(seednum), 
         yesno = as.integer(seednum>0)) %>% # was a capsule formed?
  add_combos()
  
# germination -------------------------------------------------------------

germination <- read_sheet(gsheet, "germination", col_types="c") %>% 
  mutate(across(c(planted, germ), as.integer)) %>% 
  group_by(crossid, mompop, momid, species, dadpop, dadid, momsp, crosstype) %>% 
  summarize(planted = sum(planted), germ = sum(germ), .groups="drop") %>% # add together all the pots of one cross
  mutate(vp = germ/planted) %>% # proportion of planted seeds that germinated (viable)
  add_combos()
  

