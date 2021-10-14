library(tidyverse)
library(lubridate)
library(googlesheets4)

# common ------------------------------------------------------------------

gsheet <- gs4_get("1pYbAnEDw2KfM34l85wlJV6pfAr1DroPj_7GjfApnCq8")
sheet_names(gsheet)

crosscol <- c("green","blue","orange","red")
hk.pops <- set_names(c("WK","WK","WK","879WKG","892WKG","904WPG","3587WP"),
                     c("794","866","899","879","892","904","3587")) #merge Waianae Kai populations (WK)
hk.species <- set_names(c(rep("hookeri",4),rep("kaalae",3)), names(hk.pops))

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

# crosses -----------------------------------------------------------------

crosses <- read_sheet(gsheet, "crosses", col_types="c") %>% 
  mutate(species=recode(momsp, kaal="kaalae", hook="hookeri")) %>% 
  add_combos()

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

# vegbiomass --------------------------------------------------------------

first_planting <- ymd("2016-03-10") #verified this is the first seed planting date for this set

vegbiomass <- read_sheet(gsheet, "vegbiomass", col_types="c") %>% 
  drop_na(crossid) %>% 
  mutate(across(contains("date"), ymd),
         #TODO from hybridbiomass.Rmd - is this the median planting date?
         collect = collect.date - first_planting, 
         mass=as.numeric(mass)) %>%
  left_join(crosses) %>% 
  add_combos()

# survflr -------------------------------------------------------------

survflr <- read_sheet(gsheet, "survflr", col_types="c") %>% 
  filter(crossid != "107", plantid != "0") %>% # 107 not listed in crosses, plantid=0: seeds did not germinate
  mutate(across(matches("date"), ymd),
         across(all_of(c("delay","firstinflo.biomass.mg")), as.numeric),
         alive =    ifelse(use.alive.flowered == "yes", alive == "yes", NA), #change usable statuses to boolean
         flowered = ifelse(use.alive.flowered == "yes", na_if(flowered, "?") == "yes", NA), 
         firstinflo.biomass.mg = ifelse(use.fib == "yes", firstinflo.biomass.mg, NA),
         firstflower.date = if_else(use.firstflower == "yes", firstflower.date, as.Date(NA)),
         firstflower = firstflower.date - first_planting) %>% 
  left_join(crosses) %>% 
  add_combos()
#for flowering analysis, add filter(alive) %>% drop_na(flowered)

# pollen ------------------------------------------------------------------

pollen <- read_sheet(gsheet, "pollen", col_types="c") %>% 
  filter(crosstype != "field") %>% 
  mutate(across(starts_with(c("V","I"), ignore.case=F), as.integer)) %>% rowwise() %>% 
  mutate(v = sum(c_across(num_range("V", 1:10)))/10, #average counts
         i = sum(c_across(num_range("I", 1:10)))/10) %>% ungroup() %>% 
  mutate(vpf = (v * 444.4/20) * 10, #number viable per anther (5 flrs x 4 anthers)
         ipf = (i * 444.4/20) * 10, #number inviable per anther (5 flrs x 4 anthers)
         vp = vpf / (vpf+ipf), #proportion viable
         tpf = vpf + ipf, #total pollen grains per flower
         date = ymd(date)) %>% 
  separate(fullcross, sep=" x ", into=c("mompid","dadpid"), remove=F) %>% 
  separate(mompid, into=c("mompop", "momid"), extra="merge") %>% 
  separate(dadpid, into=c("dadpop", "dadid"), extra="merge") %>% 
  mutate(species = recode(mompop, !!!hk.species),
         dadsp =   recode(dadpop, !!!hk.species),
         cross = toupper(paste0(str_sub(species,0,1), str_sub(dadsp,0,1)))) %>% 
  add_combos()

# inflobiomass ------------------------------------------------------------

inflobiomass <- read_sheet(gsheet, "inflobiomass", col_types="c") %>% 
  mutate(across(matches("date"), ymd),
         across(c(flrs, inflo.e, inflo, mass), as.numeric),
         collect = collect.date - first_planting) %>% 
  left_join(crosses) %>% 
  add_combos()

inflobiomass.sum <- inflobiomass %>% #add together envelopes for regression and the rest
  group_by(across(all_of(c("plantid", colnames(crosses))))) %>% 
  summarize(inflo = sum(inflo, na.rm=T),
            mass =  sum(mass,  na.rm=T),
            collect = max(collect, na.rm=T), .groups="drop") 

# f1seeds -----------------------------------------------------------------


# f1seedmass --------------------------------------------------------------


# f2pollen ----------------------------------------------------------------


# f2seeds -----------------------------------------------------------------


