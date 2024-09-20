source(here::here("code/01_data-cleaning.R"))
i_am("code/02_data-analysis.R")

## plot assimilation effiency by stoichio

assimilation %>%
  filter(!is.na(stoich)) %>%
  ggplot()+
  geom_line(aes(x = log10(stoich), y = assimilation_eff,  color = stoich_var, group = sourceID),
            stat = 'smooth', method = 'lm', se = FALSE, alpha = 0.5, linewidth = 1)+
  geom_point(aes(x = log10(stoich), y = assimilation_eff, group = sourceID, shape = sourceID, fill = stoich_var),
             color = 'black',size = 2)+

  guides(color = 'none',
         fill = 'none')+
  scale_shape_manual(values = c(21:25))+
  scale_y_continuous(name = "Assimilation Efficiency (%)", limits = c(0,1), expand = c(0,0.05))+
  scale_x_continuous(name = "C:X (molar)")+
  facet_wrap(~stoich_var, scales = 'free')

# growth efficiency

growth_eff %>%
  filter(!is.na(stoich)) %>%
  ggplot()+
  geom_line(aes(x = stoich, y = eff,  color = eff_var, group = sourceID),
            stat='smooth',method = 'lm',se = FALSE, alpha = 0.5,
            linewidth = 1)+
  geom_point(aes(x = stoich, y = eff, group = sourceID, shape = eff_var, fill = eff_var),
             color = 'black',size = 2)+
  # guides(color = 'none')+
  scale_shape_manual(values = 21:23)+
  scale_y_continuous(name = "Growth Efficiency (%)", limits = c(0,NA), expand = c(0,0.05))+
  scale_x_continuous(name = "C:P (molar)", limits = c(0,NA))+
  facet_wrap(~eff_var, scales = 'free')

#standardize the effects by growth variable

growth %>%
  # group_by(growth_var) %>%
  # mutate(growth_s = scale(growth, center = TRUE, scale = TRUE)) %>%
  # filter(growth_s < 2.2) %>%
  mutate(growth_s = growth) %>%
  filter(growth_var == 'd-1') %>%
  ggplot()+
  # geom_hline(aes(yintercept = 0), color = 'darkgrey', linewidth = 1)+
  geom_line(aes(x = stoich, y = growth_s, color = stoich_var, group = sourceID),
            stat = 'smooth',se = FALSE, span = 0.999, alpha = 0.5, method = 'lm',
            linewidth = 1)+
  geom_point(aes(x = stoich, y = growth_s, shape = sourceID, fill = stoich_var),
             color = 'black',size = 2)+
  guides(fill = 'none',
         color = 'none')+
  scale_shape_manual(values = c(21:25))+
  scale_y_continuous(name = "Growth rate (d-1)", expand = c(0,0.01))+
  scale_x_continuous(name = "C:X (molar)")+
  facet_wrap(~stoich_var, scales = 'free')

# standardize the effects by ingestion variable

ingestion %>%
  filter(ingestion_var %in% c("ug_C_d-1",
                              "ug_C_ug_C_d-1",
                              "mg_mg-1_h-1",
                              "mg_C_mg-1_d-1",
                              "mg_C_mg_ind-1_d-1"),
         stoich_var != 'NPmol') %>%
  group_by(ingestion_var) %>%
  mutate(ingestion = ifelse(ingestion_var == 'mg_mg-1_h-1', ingestion/1000*24,ingestion)) %>%
  mutate(ingestion_s = ingestion) %>% #scale(ingestion, center = TRUE, scale = TRUE)) %>%
  ggplot()+
  geom_hline(aes(yintercept = 0), color = 'grey', linewidth = 1)+
  geom_line(aes(x = stoich, y = log10(ingestion_s), color = stoich_var, group = sourceID),
            stat = 'smooth', method = 'lm', se = FALSE, alpha = 0.8,
            linewidth = 1)+
  geom_point(aes(x = stoich, y = log10(ingestion_s), shape = sourceID, fill = stoich_var),
             color = 'black',size = 2)+
  guides(color = 'none',
         fill = 'none')+
  scale_shape_manual(values =21:25)+
  scale_y_continuous(name = expression(log[10]~"[Ingestion (mass * "*mass^-1*"*"*time^-1*")]"))+
  facet_wrap(~stoich_var, scales = 'free')


  ingestion %>%
    filter(ingestion_var %in% c("no_squares_hr-1_ind-1",
                                "nmol_C_ind-1_h-1",
                                "ugC_ind-1_d-1",
                                "mL_h-1_ind-1"),
                                stoich_var != 'NPmol') %>%
             group_by(ingestion_var) %>%
             mutate(ingestion_s = ingestion) %>% #scale(ingestion, center = TRUE, scale = TRUE)) %>%
             ggplot()+
             # geom_hline(aes(yintercept = 0), color = 'grey', linewidth = 1)+
             geom_line(aes(x = stoich, y = log10(ingestion_s), color = stoich_var, group = sourceID),
                       stat = 'smooth', method = 'lm', se = FALSE, alpha = 0.8,
                       linewidth = 1)+
             geom_point(aes(x = stoich, y = log10(ingestion_s), shape = sourceID, fill = stoich_var),
                        color = 'black',size = 2)+
             guides(color = 'none',
                    fill = 'none')+
             scale_shape_manual(values = 21:25)+
             scale_y_continuous(name = expression(log[10]~"[Ingestion (mass * "*ind^-1*"*"*time^-1*")]"))+
             facet_wrap(~stoich_var, scales = 'free')
