library(tidyverse)
library(lubridate)
library(googlesheets4)

# common ------------------------------------------------------------------

gsheet <- gs4_get("1pYbAnEDw2KfM34l85wlJV6pfAr1DroPj_7GjfApnCq8")
datanames <- sheet_names(gsheet)[-1]

species.pops <- tibble(
  collection = c("794","866","899","879","892","904","3587"),
  pop = c("WK","WK","WK","879WKG","892WKG","904WPG","3587WP"), #merge Waianae Kai populations (WK)
  sp = c(rep("hook",4),rep("kaal",3)))

hk.pops <- set_names(species.pops$pop, species.pops$collection) 
hk.species <- set_names(rep(species.pops$sp, 2), c(species.pops$collection, species.pops$pop))

first_planting <- ymd("2016-03-10") #verified this is the first seed planting date for this set

sxc.levels <- c("kaal.within", "kaal.between", "kaal.hybrid", "hook.hybrid", "hook.between", "hook.within")

pop_factors <- function(data) {
  mutate(data,
         across(where(is.character), as.factor),
         across(ends_with("pop"), recode_factor, !!!hk.pops),
         across(any_of("sxc"), fct_relevel, sxc.levels))
}

add_combos <- function(data) {
  mutate(data, 
         sxc = paste(momsp, crosstype, sep="."),
         sxcxm = paste(momsp, crosstype, mompop, momid, sep="."),
         mompid = paste(mompop, momid, sep="."),
         across(any_of("dadid"), ~paste(dadpop, ., sep="."), .names="dadpid"),
         smompop = paste(momsp, mompop, sep="")) %>% 
    pop_factors() %>% 
    mutate(crosstype = fct_relevel(crosstype, "hybrid", after = 2))
}

split_full_crosses <- function(data.full) {
  data.mom <- bind_rows(
    .id = "momgeneration", 
    F = data.full %>% filter(str_detect(momfullcross,"x")) %>% 
      separate(momfullcross, sep=" x ", into=c("mommompid","momdadpid"), remove=F) %>% 
      separate(mommompid, into=c("mommompop", "mommomid"), extra="merge") %>% 
      separate(momdadpid, into=c("momdadpop", "momdadid"), extra="merge"),
    P = data.full %>% filter(!str_detect(momfullcross,"x")) %>% 
      separate(momfullcross, into=c("mompop", "momid"), extra="merge", remove=F))
  
  data.mom.dad <- list(
    F = data.mom %>% filter(str_detect(dadfullcross,"x")) %>% 
      separate(dadfullcross, sep=" x ", into=c("dadmompid","daddadpid"), remove=F) %>% 
      separate(dadmompid, into=c("dadmompop", "dadmomid"), extra="merge") %>% 
      separate(daddadpid, into=c("daddadpop", "dadadid"), extra="merge"))
  
  if(("crosstype" %in% colnames(data.full)) && ("control" %in% data.full$crosstype)) {
    data.mom.dad[["control"]] <- data.mom %>% filter(crosstype == "control") %>% mutate(dadfullcross=NA)
    data.mom.dad[["P"]] <- data.mom %>% filter(!str_detect(dadfullcross,"x"), crosstype != "control")
  } else {
    data.mom.dad[["P"]] <- data.mom %>% filter(!str_detect(dadfullcross,"x"))
  }
  data.mom.dad[["P"]] <- data.mom.dad[["P"]] %>% separate(dadfullcross, into=c("dadpop", "dadid"), extra="merge", fill="right", remove=F)
  
  bind_rows(data.mom.dad, .id = "dadgeneration") %>% 
    mutate(generation = paste0(momgeneration,dadgeneration))
}

initials <- function(a, b="") toupper(paste0(str_sub(a,0,1), str_sub(b,0,1)))

pop_to_species <- function(data) {
  mutate(data, 
         across(ends_with("pop"), ~ recode(.x, !!!hk.species), 
                .names="{str_remove(.col,\"pop\")}sp"),
         momcross = ifelse(is.na(momsp), initials(mommomsp, momdadsp), initials(momsp)) %>% 
           factor(levels=c("H",  "HK", "KH", "K")),
         dadcross = ifelse(is.na(dadsp), initials(dadmomsp, daddadsp), initials(dadsp)) %>% 
           recode(NANA="control") %>% factor(levels=c("H",  "HK", "KH", "K", "control")) %>% droplevels(),
         momdadcross = paste(momcross, dadcross, sep = " x "))
}

calc_pollen <- function(data) {
  mutate(data, across(starts_with(c("V","I"), ignore.case=F), as.integer)) %>% rowwise() %>% 
    mutate(v = sum(c_across(num_range("V", 1:10)))/10, #average counts
           i = sum(c_across(num_range("I", 1:10)))/10) %>% ungroup() %>% 
    mutate(vpf = (v * 444.4/20) * 10, #number viable per anther (5 flrs x 4 anthers)
           ipf = (i * 444.4/20) * 10, #number inviable per anther (5 flrs x 4 anthers)
           vp = vpf / (vpf+ipf), #proportion viable
           tpf = vpf + ipf, #total pollen grains per flower
           date = ymd(date))
}

# crosses -----------------------------------------------------------------

crosses <- read_sheet(gsheet, "crosses", col_types="c") %>% 
  add_combos()

# seeds -------------------------------------------------------------------

seeds <- read_sheet(gsheet, "seeds", col_types="c") %>% 
  filter(dadpop != "closed") %>% 
  mutate(viable.seeds = as.integer(viable.seeds), 
         capsule.formed = as.integer(viable.seeds>0)) %>% 
  add_combos()
  
# germination -------------------------------------------------------------

germination <- read_sheet(gsheet, "germination", col_types="c") %>% 
  mutate(across(c(planted, germinated), as.integer)) %>% 
  group_by(crossid, mompop, momid, momsp, dadpop, dadid, dadsp, crosstype) %>% 
  summarize(planted = sum(planted), germinated = sum(germinated), .groups="drop") %>% # add together all the pots of one cross
  mutate(prop.germ = germinated/planted) %>% # proportion of planted seeds that germinated 
  add_combos()

# vegbiomass --------------------------------------------------------------

vegbiomass <- read_sheet(gsheet, "vegbiomass", col_types="c") %>% 
  drop_na(crossid) %>% 
  mutate(across(contains("date"), ymd),
         collect = collect.date - first_planting, 
         veg.biomass.g=as.numeric(veg.biomass.g)) %>%
  left_join(crosses) %>% 
  drop_na(crosstype) %>% 
  add_combos()

# survflr -------------------------------------------------------------

survflr <- read_sheet(gsheet, "survflr", col_types="c") %>% 
  filter(crossid != "107", plantid != "0") %>% # 107 not listed in crosses, plantid=0: seeds did not germinate
  mutate(across(matches("date"), ymd),
         across(all_of(c("delay","firstinflo.biomass.mg")), as.numeric),
         firstinflo.biomass.g = firstinflo.biomass.mg/1000,
         alive =    ifelse(use.alive.flowered == "yes", alive == "yes", NA), #change usable statuses to boolean
         flowered = ifelse(use.alive.flowered == "yes", na_if(flowered, "?") == "yes", NA), 
         firstinflo.biomass.g = ifelse(use.fib == "yes", firstinflo.biomass.g, NA),
         collect = firstinflo.collect.date - first_planting, 
         firstflower.date = if_else(use.firstflower == "yes", firstflower.date, as.Date(NA)),
         firstflower = firstflower.date - first_planting) %>% select(-firstinflo.biomass.mg) %>% 
  left_join(crosses) %>% 
  add_combos()
#for flowering analysis, add filter(alive) %>% drop_na(flowered)

# pollen ------------------------------------------------------------------

pollen <- read_sheet(gsheet, "pollen", col_types="c") %>% 
  filter(crosstype != "field") %>% 
  calc_pollen() %>% 
  separate(fullcross, sep=" x ", into=c("mompid","dadpid"), remove=F) %>% 
  separate(mompid, into=c("mompop", "momid"), extra="merge") %>% 
  separate(dadpid, into=c("dadpop", "dadid"), extra="merge") %>% 
  mutate(momsp = recode(mompop, !!!hk.species),
         dadsp = recode(dadpop, !!!hk.species),
         cross = toupper(paste0(str_sub(momsp,0,1), str_sub(dadsp,0,1)))) %>% 
  add_combos()

#TODO investigate if crossid in id column matches up with fullcross

# inflobiomass ------------------------------------------------------------

inflobiomass <- read_sheet(gsheet, "inflobiomass", col_types="c") %>% 
  mutate(across(matches("date"), ymd),
         across(c(flrs, inflo.e, inflo, inflo.biomass.g), as.numeric),
         collect = collect.date - first_planting) %>% 
  left_join(crosses) %>% 
  add_combos()

inflobiomass.sum <- inflobiomass %>% #add together envelopes for regression and the rest
  group_by(across(all_of(c("plantid", colnames(crosses))))) %>% 
  summarize(inflo = sum(inflo, na.rm=T),
            inflo.biomass.g =  sum(inflo.biomass.g,  na.rm=T),
            collect = suppressWarnings(max(collect, na.rm=T)) %>% na_if(-Inf), .groups="drop") 

# f1seeds -----------------------------------------------------------------

f1seeds <- read_sheet(gsheet, "f1seeds", col_types="c") %>% 
  filter(dadfullcross != "closed") %>% #not sure what closed means
  mutate(viable.seeds=as.integer(viable.seeds)) %>% 
  split_full_crosses() %>% 
  pop_to_species() %>% 
  pop_factors()

# f1seedmass --------------------------------------------------------------

f1seedmass <- read_sheet(gsheet, "f1seedmass", col_types="c") %>% 
  mutate(tenseeds.mass.mg=as.numeric(tenseeds.mass.mg), #mass of 10 seeds
         seed.mass.mg = tenseeds.mass.mg / 10) %>% 
  split_full_crosses() %>% 
  pop_to_species() %>% 
  pop_factors()

# f2pollen ----------------------------------------------------------------

f2pollen <- read_sheet(gsheet, "f2pollen", col_types="c") %>% 
  calc_pollen() %>% 
  split_full_crosses() %>% 
  pop_to_species() %>% 
  pop_factors()

# f2seeds -----------------------------------------------------------------

f2seeds <- read_sheet(gsheet, "f2seeds", col_types="c") %>% 
  mutate(across(contains("date"), ymd),
         across(contains("seeds"), as.integer),
         capsule.formed = viable.seeds > 0) %>% 
  split_full_crosses() %>% 
  pop_to_species() %>% 
  pop_factors()

# export ------------------------------------------------------------------
walk(datanames, ~write_tsv(get(.), paste0("data/",., ".tsv")))

save.image("data/rims_data.rda")



