# Teste do Modelo de Jogos
# Equilíbrio de Mercado

# Pacotes
if(!require(tidyverse)){install.packages("tidyverse");require(tidyverse)}
if(!require(nloptr)){install.packages("nloptr");require(nloptr)}

# Detalhamentos:
# Todos os dados do Banco Mundial são do ano de 2021. 
# As distâncias são do CEPII

# Number of countries:
N <- 218

# Manipulação
In=diag(N)
i_n=matrix(rep(1,N),ncol=1)

# W2: Impacto homogêneo do desmatamento entre os países
W2 <- matrix(rep(1,N*N),ncol=N)
# Suponha agora que o impacto de W2 sobre a utilidade dependa positivamente do nível de florestas
forest_2021 <- WDI::WDI(indicator = "AG.LND.FRST.K2" ,start=2021,end=2021)
country_codes <- read_csv("https://gist.githubusercontent.com/tadast/8827699/raw/61b2107766d6fd51e2bd02d9f78f6be081340efc/countries_codes_and_coordinates.csv") %>% 
  select(`Alpha-2 code`,`Alpha-3 code`) %>% group_by(`Alpha-3 code`) %>% 
  summarise(`Alpha-2 code`=first(`Alpha-2 code`)) %>% ungroup()
country_sample <- readxl::read_xls("data/dist_cepii.xls") %>% 
  select('iso_o') %>% unique() %>% arrange()
country_iso2 <- country_sample %>% inner_join(country_codes,by=c("iso_o"='Alpha-3 code')) %>% na.omit()

# Agora, podemos plugar esses valores na matriz de florestas:
forests_selected <- country_iso2 %>% select(`Alpha-2 code`) %>% 
  left_join(forest_2021,by=c(`Alpha-2 code`='iso2c')) %>% 
  select(`AG.LND.FRST.K2`) %>% 
  rename(forests=`AG.LND.FRST.K2`) %>% 
  mutate(forests=if_else(is.na(forests),0,forests),
         forests=forests/sum(forests)) %>% as.matrix()
# Criando a matrix agora:
W2 <- matrix(rep(forests_selected, each = N), nrow = nrow(forests_selected))

gdp <- WDI::WDI(indicator = "NY.GDP.PCAP.KD" ,start=2021,end=2021)
gdp_selected <- country_iso2 %>% select(`Alpha-2 code`) %>% 
  left_join(gdp,by=c(`Alpha-2 code`='iso2c')) %>% 
  select(`NY.GDP.PCAP.KD`) %>% 
  rename(gdp=`NY.GDP.PCAP.KD`) %>% 
  mutate(gdp=if_else(is.na(gdp),0,gdp)) %>% as.matrix()

# NY.GDP.PCAP.KD: GDP per capita (constant 2015 US$)

# W1: Distâncias entre os países:
dist_inv <- readxl::read_xls("data/dist_cepii.xls") %>% 
  filter(iso_o %in% country_iso2$iso_o & iso_d %in% country_iso2$iso_o) %>% 
  select(iso_o,iso_d,dist) %>% pivot_wider(names_from = 'iso_d',values_from = 'dist') %>% 
  select(-iso_o) %>% 
  mutate_all(~ ifelse(. == 0, 0, 1/.)) %>% 
  as.matrix()
W1 <- dist_inv%*%(gdp_selected%*%t(gdp_selected))
# Pensando em ponderar a distância pelos PIBs dos países num esquema bem gravitacional

# Modelo de equilíbrio de mercado
eq1_compet <- function(tau){
  return(sum((In-W1%*%kronecker(tau,t(i_n))%*%In-sweep(W1, 2, tau, "*"))^2))
}
# resolvendo o modelo
x0 <- runif(N)
opts <- list("algorithm"="NLOPT_LN_SBPLX",
             "xtol_rel"=1.0e-8)
res <- nloptr( x0=x0, 
               eval_f=eq1_compet,
               lb=rep(0,N),
               opts=opts)
sol_market_tau <- res$solution

# Desmatamento
eq2_compet <- function(d){
  return(sum((In-W2%*%kronecker(d,t(i_n))%*%In-sweep(W2, 2, d, "*"))^2))
}
# resolvendo o modelo
x0 <- runif(N)
opts <- list("algorithm"="NLOPT_LN_NELDERMEAD",
             "xtol_rel"=1.0e-8)
res <- nloptr( x0=x0, 
               eval_f=eq2_compet,
               lb=rep(0,N),
               opts=opts)
sol_market_d <- res$solution

# Planejador Central
eq1_central <- function(tau){
  sum <- rep(0,N)
  db <- rep(0,N)
  da <-rep(1,N)
  for (i in 1:N){
    sum[i] <- sum(W1[,i])
    db[i] <- sum(tau*sum[i])
  }
  return(sum((2*db-da)^2))
}
# Outro solver
require(nloptr)
x0 <- runif(N)+sol_market_tau
opts <- list("algorithm"="NLOPT_LN_NELDERMEAD",
             "xtol_rel"=1.0e-8)
res <- nloptr( x0=x0, 
               eval_f=eq1_central,
               lb=rep(0,N),
               opts=opts)
sol_central_tau <- res$solution

# Desmatamento
eq2_central <- function(d){
  sum <- rep(0,N)
  db <- rep(0,N)
  da <-rep(1,N)
  for (i in 1:N){
    sum[i] <- sum(W2[,i])
    db[i] <- sum(d*sum[i])
  }
  return(sum((2*db-da)^2))
}
# Outro solver
require(nloptr)
x0 <- pmin(runif(N),sol_market_d)
opts <- list("algorithm"="NLOPT_LN_NELDERMEAD",
             "xtol_rel"=1.0e-8)
res <- nloptr( x0=x0, 
               eval_f=eq2_central,
               lb=rep(0,N),
               opts=opts)
sol_central_d <- res$solution

# Comparação entre soluções de mercado e central
sol <- cbind(sol_market_tau,sol_central_tau,sol_market_d,sol_central_d) %>% 
  as_tibble() %>% 
  mutate(tau_maior_central=sol_market_tau>=sol_central_tau,
         d_maior_mercado=sol_market_d>=sol_central_d)

# Analisando os resultados dos países com mais florestas no modelo:
forests_vec <- c(168,28,34,206,38)
# Russia, Brazil, Canada, United States of America e China
sol[forests_vec,]

# Outras matrizes W1
#### CORRENTE COMERCIAL ##### 
corrente_comercial <- WDI::WDI(indicator = "NE.TRD.GNFS.ZS" ,start=2021,end=2021)
corrente_comercial_selected <- country_iso2 %>% select(`Alpha-2 code`) %>% 
  left_join(corrente_comercial,by=c(`Alpha-2 code`='iso2c')) %>% 
  select(`NE.TRD.GNFS.ZS`) %>% 
  rename(trade_share=`NE.TRD.GNFS.ZS`) %>% 
  mutate(trade_share=if_else(is.na(trade_share),0,trade_share)) %>% as.matrix()
W1 <- matrix(rep(corrente_comercial_selected,each=N),nrow = nrow(corrente_comercial_selected))

#### ÍNDICE DE COMPLEXIDADE ECONOMICA #####
indice_complexidade_economica <- read_csv("data/indice_complexidade_economica.csv") %>% 
  select(Country,`ECI 2021`)
country_codes1 <- read_csv("https://gist.githubusercontent.com/tadast/8827699/raw/61b2107766d6fd51e2bd02d9f78f6be081340efc/countries_codes_and_coordinates.csv")
gdp <- WDI::WDI(indicator = "NY.GDP.PCAP.KD" ,start=2021,end=2021)

country_code2 <- country_codes1 %>% select(Country,`Alpha-2 code`)
df1 <- gdp %>% select(`iso2c`,`NY.GDP.PCAP.KD`) %>% 
  inner_join(country_code2,by=c('iso2c'="Alpha-2 code")) %>% 
  fuzzyjoin::stringdist_inner_join(indice_complexidade_economica,by=c("Country"),
                                   method="lv",
                                   max_dist=2) %>% 
  select(-starts_with("Country")) %>% 
  left_join(country_codes3,by=c('iso2c'="Alpha-2 code")) %>% 
  group_by(`Alpha-3 code`) %>% 
  summarise(iso2c=first(iso2c),
            `NY.GDP.PCAP.KD`=first(`NY.GDP.PCAP.KD`),
            `ECI 2021`=first(`ECI 2021`))
country_codes3 <- country_codes1 %>% select(`Alpha-2 code`,`Alpha-3 code`)

dist_inv1 <- readxl::read_xls("data/dist_cepii.xls") %>% 
  filter(iso_o %in% df1$`Alpha-3 code` & iso_d %in% df1$`Alpha-3 code`) %>% 
  select(iso_o,iso_d,dist) %>% pivot_wider(names_from = 'iso_d',values_from = 'dist') 
dist_inv1_codes <- dist_inv1$iso_o
dist_inv12 <- dist_inv1 %>% 
  select(-iso_o) %>% 
  mutate_all(~ ifelse(. == 0, 0, 1/.)) %>% 
  as.matrix()

df_codes_selected <- df1 %>% inner_join(dist_inv1,by=c("Alpha-3 code"='iso_o')) %>% 
  select(`iso2c`)
forests_selected1 <- df_codes_selected %>% 
  left_join(forest_2021,by=c('iso2c')) %>% 
  select(`AG.LND.FRST.K2`) %>% 
  rename(forests=`AG.LND.FRST.K2`) %>% 
  mutate(forests=if_else(is.na(forests),0,forests),
         forests=forests/sum(forests)) %>% 
  as.matrix()

# GDP, Florestas, Distâncias e ECI calculados
gdp_selected <- df1 %>% filter(iso2c %in% df_codes_selected$iso2c) %>% 
  select(`NY.GDP.PCAP.KD`)
eci_2021 <- df1 %>% filter(iso2c %in% df_codes_selected$iso2c) %>% 
  select(`ECI 2021`)
eci_2021 <- matrix(eci_2021$`ECI 2021`,ncol=1)
dist_inv2 <- as.matrix(dist_inv1[,2:132])
W1 <- dist_inv2%*%(eci_2021%*%t(eci_2021))
W2 <- matrix(rep(forests_selected1, each = nrow(forests_selected1)), nrow = nrow(forests_selected1))

# Resultados
forests_vec <- c(104,18,20,23)
sol[forests_vec,]

# Salvando as matrizes de dados
writexl::write_xlsx(as.data.frame(W1),"W1_gdp.xlsx")
writexl::write_xlsx(as.data.frame(W2),"W2.xlsx")
writexl::write_xlsx(as.data.frame(W1),"W1_trade.xlsx")
writexl::write_xlsx(as.data.frame(W1),"W1_eic.xlsx")
writexl::write_xlsx(as.data.frame(W2),"W2_eic.xlsx")

# Salvando os resultados
writexl::write_xlsx(as.data.frame(sol),"sol_gdp_percapita.xlsx")
writexl::write_xlsx(as.data.frame(sol),"sol_corrente_comercial.xlsx")
writexl::write_xlsx(as.data.frame(country_iso2),"ordem_paises.xlsx")
