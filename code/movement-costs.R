library(tidyverse)
theme_set(theme_minimal())

Mb = seq(0.0001, 10000, length.out = 100)
v = seq(0.001, 2, length.out = 100)
e_Mb = 10.7*(Mb^-0.316)*v + 6.03*Mb^(-0.303)

#xpand grid to view patterns
pred_grid = expand.grid(Mb = Mb,v = v) %>%data.frame %>%  mutate(e_Mb = 10.7*Mb^-0.316*v + 6.03*Mb^-0.303)

pred_grid %>%
  ggplot()+
  geom_line(aes(x = Mb, y = e_Mb, group = v, color = v))+
  viridis::scale_color_viridis()+
  scale_y_log10(name = expression(log[10]*"( "*frac("E",M[b])*" )"~"[W "*kg^-1*"]"))+
  scale_x_log10(name = expression(log[10]*"( "*M[b]*" ) [kg]"))+
  labs(color = expression("v [ m"~s^-1~"]"))

# total energy to move 1 m based on a bunch of scaling
pred_grid = pred_grid %>% mutate(eJ_m = 220*Mb^0.485,
                                 eJ_m_Mb = eJ_m/Mb,
                                 eJ_m_l = eJ_m/(Mb^0.66))

pred_grid %>%
  ggplot()+
  geom_line(aes(x = Mb, y = eJ_m))+
  scale_y_continuous(name = expression("Energy cost of travel (J "~kg^-1*" )"))+
  scale_x_continuous(name = "Body size (kg)")


pred_grid %>%
  ggplot()+
  geom_line(aes(x = Mb, y = eJ_m_l))+
  scale_y_log10(name = expression("Energy cost of travel (J "~kg^-1*" )"))+
  scale_x_log10(name = expression("Body surface area ( "*m^2*" )"))
