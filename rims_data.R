library(tidyverse)
library(googlesheets4)

# common ------------------------------------------------------------------

gs4_auth(email = T)
traits <- read_sheet("1f487VMTsYzYGn0_daogtTF0nNE7ntoqd0Yh4h9tgyDE", "Traits") #RIMs summary

gsheet <- gs4_get("1pYbAnEDw2KfM34l85wlJV6pfAr1DroPj_7GjfApnCq8") #RIMs data
datanames <- sheet_names(gsheet)[-c(1,2,4)]#all except crosses, seeds_full, germination_full

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
    mutate(viable = sum(c_across(num_range("V", 1:10))), #total counts
           inviable = sum(c_across(num_range("I", 1:10)))) %>% ungroup() %>% 
    mutate(total = viable + inviable,
           across(c(viable, inviable, total), ~ . * 444.4 / (5 * 4), .names="{col}.per.anther"), #444.4 from hemocytometer, 5 flrs, 4 anthers
           prop.viable = viable / (viable + inviable), 
           date = ymd(date))
}

# crosses -----------------------------------------------------------------

crosses <- read_sheet(gsheet, "crosses", col_types="c") %>% 
  add_combos()

# seeds -------------------------------------------------------------------

seeds_full <- read_sheet(gsheet, "seeds_full", col_types="c", na=c("NA","")) %>% 
  mutate(across(ends_with("date"), ymd)) %>% 
  separate(mompop.type, c("momsp","momlocality"), sep=" ") %>% 
  separate(dadpop.type, c("dadsp","dadlocality"), sep=" ")
  
seeds <- read_sheet(gsheet, "seeds", col_types="c") %>% 
  filter(dadpop != "closed") %>% #not sure what closed means
  mutate(viable.seeds = as.integer(viable.seeds), 
         capsule.formed = as.integer(viable.seeds>0)) %>% 
  left_join(select(seeds_full, labelid, dadid)) %>% 
  add_combos()

seeds.nonzero <- filter(seeds, viable.seeds > 0)
  
# germination -------------------------------------------------------------

germination <- read_sheet(gsheet, "germination") %>% 
  group_by(crossid, mompop, momid, momsp, dadpop, dadid, dadsp, crosstype) %>% #drop potid, comments
  summarize(planted = sum(planted), germinated = sum(germinated), .groups="drop") %>% # add together all the pots of one cross
  mutate(prop.germ = germinated/planted) %>% # proportion of planted seeds that germinated 
  add_combos() %>% 
  left_join(#estimate date that half of seedlings have germinated
    read_sheet(gsheet, "germination_full") %>% 
      filter(crosstype!="control") %>% #one pot with seeds from a control (unpollinated) cross, not present in "germination" sheet
      pivot_longer(starts_with("2016"), names_to = "scoredate", values_to = "germinated") %>% 
      group_by(crossid, tray, plantdate, scoredate) %>% #drop pot identifiers, comments.germ, and ID columns updated in "germination" sheet
      summarize(across(c(planted, germinated), sum), .groups = "drop") %>% 
      mutate(scoredate = ymd(scoredate), scoreday = as.numeric(scoredate - ymd(as.character(plantdate)), units="days")) %>% 
      group_by(crossid) %>% mutate(prop.germ = germinated/max(germinated)) %>% #scale to maximum seedlings, not seeds planted
      filter(!all(germinated == 0), #exclude crosses that never germinated 
             !all(prop.germ == 1)) %>%  #exclude crosses that had already all germinated
      mutate(scoreday.peak = scoreday[min(which(prop.germ==1))]) %>% #day that seedling number is maximum
      mutate(prop.germ = if_else(scoreday > scoreday.peak, 1, prop.germ)) %>% # erase death after peak to focus on timing
      ungroup() %>% drop_na(prop.germ) %>% nest(.by=crossid) %>% 
      mutate(coef = map(data, ~coef(glm(prop.germ~scoreday, family="quasibinomial", data=.x))), .keep="unused") %>% 
      unnest_wider(coef) %>% mutate(half.germ.day = -`(Intercept)`/scoreday, .keep="unused"))

# vegbiomass --------------------------------------------------------------

vegbiomass <- read_sheet(gsheet, "vegbiomass", col_types="c") %>% 
  filter(is.na(use.veg)) %>%  
  mutate(across(contains("date"), ymd),
         collect = collect.date - first_planting, 
         veg.biomass.g=as.numeric(veg.biomass.g)) %>%
  left_join(crosses) %>% 
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
         firstflower = as.integer(firstflower.date - first_planting)) %>% select(-firstinflo.biomass.mg) %>% 
  left_join(crosses) %>% 
  add_combos()

survflr.alive <- survflr %>% filter(alive)

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
  filter(is.na(use.inflo)) %>% 
  mutate(across(matches("date"), ymd),
         across(c(flrs, inflo.e, inflo, inflo.biomass.g), as.numeric),
         collect = as.numeric(collect.date - first_planting, units="days"),
         flrs.per.inflo.biomass.g = flrs/inflo.biomass.g) %>% 
  left_join(crosses) %>% 
  add_combos()

biomass.flrs.lm <- lm(log10(flrs) ~ poly(log10(inflo.biomass.g),2) * cross, data=inflobiomass)

inflobiomass.sum <- inflobiomass %>% #add together envelopes for regression and the rest
  group_by(across(all_of(c("plantid", colnames(crosses))))) %>% 
  summarize(inflo = sum(inflo, na.rm=T),
            tot.inflo.biomass.g =  sum(inflo.biomass.g,  na.rm=T),
            inflo.biomass.g = tot.inflo.biomass.g / inflo, #need to feed the average biomass per inflo into the model
            collect = suppressWarnings(max(collect, na.rm=T)) %>% na_if(-Inf),
            flrs = median(flrs, na.rm=T),
            flrs.per.inflo.biomass.g = median(flrs.per.inflo.biomass.g, na.rm=T),
            .groups="drop") %>% 
  mutate(flower.number = inflo * 10 ^ predict(biomass.flrs.lm, newdata = .), #multiply by inflo count *after* prediction
         avg.inflo.biomass.g = inflo.biomass.g,
         inflo.biomass.g = tot.inflo.biomass.g) #the name used in analysis for total inflo biomass

# f1seeds -----------------------------------------------------------------

f1seeds <- read_sheet(gsheet, "f1seeds", col_types="c") %>% 
  filter(dadfullcross != "closed") %>% #not sure what closed means
  mutate(viable.seeds = as.integer(viable.seeds), 
         capsule.formed = viable.seeds > 0) %>% 
  split_full_crosses() %>% 
  pop_to_species() %>% 
  filter(dadcross != "control", momdadcross!="KH x HK") %>% #single cross of that type
  pop_factors()

f1seeds.nonzero <- filter(f1seeds, viable.seeds > 0)

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
