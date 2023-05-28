# load libraries 
library(lessR)
library(cld3)
library(quanteda)
library(dplyr)
library(ggplot2)
library(stm)
library(reshape2)
library(stminsights)
library(ggthemes)
library(tidyverse)
library(ggraph)

# load data 
WA_data = read.csv("/Users/tatiana/Desktop/China_2023/pub_pics/Suppl/Data.csv",  sep=";")
str(WA_data)

# create bar plot for paper distribution over corresponding author country
BarChart(COUNTRY, data=WA_data, fill="viridis", ylab="Number of articles", sort="+")

# create bar  plot for paper distribution over year of publication
BarChart(YEAR, data=WA_data, fill=(count), xlab="Year", ylab="Number of articles")

# remove copyright from abstracts
WA_data$ABSTRACT=gsub("©.*","", WA_data$ABSTRACT)

# add new variable to dataset

WA_data$TEXT=WA_data$ABSTRACT

# remove non-English abstracts
WA_data=subset(WA_data, detect_language(WA_data$TEXT) == "en")
length(WA_data$TEXT)
dim(WA_data)

# create corpus and docvars
WA_corp1 <- corpus(WA_data, text_field = "TEXT", 
                   docid_field = "doc_id", unique_docnames = TRUE)
docvars(WA_corp1,  "text") <- WA_data$TEXT  
 

# preprocess corpus 
WA_tokens1 <- tokens(WA_corp1, remove_numbers = T, remove_punct = T, remove_symbols = T)
WA_tokens1 <- WA_tokens1 %>% 
  tokens_tolower() %>% 
  tokens_remove(stopwords('english'), padding = TRUE)
WA_tokens1

# find collocations
colls <- textstat_collocations(WA_tokens1,
                               min_count = 50)

 
#add collocations to tokens

WA_tokens1 <- tokens_compound(WA_tokens1, colls, join = TRUE) %>% 
  tokens_remove('')  
WA_tokens1 


## create Document-Feature Matrix (DFM)
dfm_WA1 <- dfm(WA_tokens1)
dim(dfm_WA1)
dfm_WA1

## inspect DFM 
tstat_freq <- textstat_frequency(dfm_WA1, n = 50)
head(tstat_freq, 100)

## lowering DFM sparsity 
dfm_WA2 <- dfm_WA1 %>% 
  dfm_keep(min_nchar = 2)  
  dfm_trim(min_docfreq = 0.01, max_docfreq = 0.50,  
           docfreq_type = 'prop')

# convert DFM  to  STM format  
out <- convert(dfm_WA2, to = 'stm')
names(out)

# choosing number of topics

K<-seq(10,50, by=5) 
kresult <- searchK(out$documents, out$vocab, K,  data=out$meta, prevalence =~ s(YEAR) + COUNTRY, verbose=T)
kresult

# inspect and plot topic modeling quality metrics  
print(kresult$results)
options(repr.plot.width=6, repr.plot.height=6)
plot(kresult)

# inspect tradeoff between coherence and exclusivity 
plot <- data.frame("K" =K, 
                   "Coherence" = unlist(kresult$results$semcoh),
                   "Exclusivity" = unlist(kresult$results$exclus))

# plot tradeoff between coherence and exclusivity

plot <- melt(plot, id=c("K"))
ggplot(plot, aes(K, value, color = variable)) +
  geom_line(size = 1.5, show.legend = FALSE) +
  facet_wrap(~variable,scales = "free_y") +
  labs(x = "Number of topics")

# build and inspect three models for final choice (K=25, 30, 35)
# 25-topic model
mod.25 <- stm::stm(out$documents, 
                  out$vocab, K=25, 
                  data=out$meta, 
                  prevalence = ~s(YEAR) + COUNTRY,
                  verbose = T)
summary(mod.25)

# 30-topic model
mod.30 <- stm::stm(out$documents, 
                   out$vocab, K=30, 
                   data=out$meta, 
                   prevalence = ~s(YEAR) + COUNTRY,
                   verbose = T)
summary(mod.30)

# 35-topic model
mod.35 <- stm::stm(out$documents, 
                   out$vocab, K=35, 
                   data=out$meta, 
                   prevalence = ~s(YEAR) + COUNTRY,
                   verbose = T)

summary(mod.35)

# Model evaluation using tradeoff method 

mod25df<-as.data.frame(cbind(c(1:25),exclusivity(mod.25), semanticCoherence(model=mod.25, out$documents), "PQ25T"))
mod30df<-as.data.frame(cbind(c(1:30),exclusivity(mod.30), semanticCoherence(model=mod.30, out$documents), "PQ30T"))
mod35df<-as.data.frame(cbind(c(1:35),exclusivity(mod.35), semanticCoherence(model=mod.35, out$documents), "PQ35T"))

models<-rbind(mod25df, mod30df, mod35df)
colnames(models)<-c("Topic","Exclusivity", "SemanticCoherence", "Model")

models$Exclusivity<-as.numeric(as.character(models$Exclusivity))
models$SemanticCoherence<-as.numeric(as.character(models$SemanticCoherence))

options(repr.plot.width=7, repr.plot.height=6, repr.plot.res=100)

plotmodels <-ggplot(models, aes(SemanticCoherence, Exclusivity, color = Model))+
  geom_point(size = 2, alpha = 0.7) + 
  geom_text(aes(label=Topic), nudge_y=.04)+
  labs(x = "Semantic coherence",
       y = "Exclusivity",
       title = "Comparing exclusivity and semantic coherence")
plotmodels

# prepare effects for YEAR
year25 <- stm::estimateEffect(1:25 ~ s(YEAR), mod.25, meta=out$meta)
year30 <- stm::estimateEffect(1:30 ~ s(YEAR), mod.30, meta=out$meta)
year35 <- stm::estimateEffect(1:20 ~ s(YEAR), mod.35, meta=out$meta)

# prepare effects for COUNTRY
country25 <- stm::estimateEffect(1:25 ~ COUNTRY, mod.25, meta=out$meta)
country30 <- stm::estimateEffect(1:30 ~ COUNTRY, mod.30, meta=out$meta)
country35 <- stm::estimateEffect(1:35 ~ COUNTRY, mod.35, meta=out$meta)

# inspect three models using library "stminsights" 
run_stminsights()


# use toLDAvis() function from the stm package to create visualizations for exploring topic and word distributions using LDAvis topic browser:
toLDAvis(mod = mod.30, docs = out$documents, reorder.topics=F)

# inspect topic distribution in document collections and words characterising topics
td_beta <- tidy(mod.30)
td_gamma <- tidy(mod.30, matrix = "gamma")

top_terms <- td_beta %>%
  arrange(beta) %>%
  group_by(topic) %>%
  top_n(7, beta) %>%
  arrange(-beta) %>%
  select(topic, term) %>%
  summarise(terms = list(term)) %>%
  mutate(terms = map(terms, paste, collapse = ", ")) %>% 
  unnest()

gamma_terms <- td_gamma %>%
  group_by(topic) %>%
  summarise(gamma = mean(gamma)) %>%
  arrange(desc(gamma)) %>%
  left_join(top_terms, by = "topic") %>%
  mutate(topic = paste0("Topic ", topic),
         topic = reorder(topic, gamma))
gamma_terms%>%print(n=30)


gamma_terms %>%
  top_n(30, gamma) %>%
  ggplot(aes(topic, gamma, label = terms, fill = topic)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 0.06),
                     labels = percent_format(suffix="")) +
  geom_text(hjust = 0, nudge_y = 0.0005, size = 3)+ 
  theme(plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 12)) +
  labs(x = NULL, y = expression(gamma),
       title = "Top 30 topics by prevalence in the Word Association Research Dataset",
       subtitle = "With the top words that contribute to each topic")

# topic distribution over years

topicprop30<-make.dt(mod.30, meta)
topic30prop <- topicprop30 %>% select(c(2:31))
topic_proportion_per_year30 <- aggregate(topic30prop, by = list(Year = WA_data$YEAR), mean)
fig30 <- melt(topic_proportion_per_year30, id.vars = "Year")

ggplot(fig30, aes(x=Year, y=value, fill=variable)) + 
  geom_bar(stat = "identity") + ylab("proportion") + 
  scale_fill_manual(values = paste0(polychrome(30), "FF"), name = "Topic", labels=labels) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Topic proportion over time in Word Association Research Dataset")

# topic distribution over countries

topic_proportion_per_country30 <- aggregate(topic30prop, by = list(Country = WA_data$COUNTRY), mean)
country30 <- melt(topic_proportion_per_country30, id.vars = "Country")

ggplot(country30, aes(x=Country, y=value, fill=variable)) + 
  geom_bar(stat = "identity") + ylab("proportion") + 
  scale_fill_manual(values = paste0(polychrome(30), "FF"), name = "Topic", labels=labels) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Topic proportion over country in Word Association Research Dataset")

#extract topic correlation network as tidygraph objects and add labels and topic proportions
stm_corrs <- get_network(model = mod.30,
                         method = 'simple',
                         cutoff=0.05,
                         labels = labels,
                         cutiso = TRUE)

labels = c("Brand research", "Text mining", "health_care", "Network_analysis", "Political_studies", "Dream_studies", "L2_studies",  "Cultural_Discourse_analysis", "Student_cognition", "Emotion_research", "Psychology", "Dictionary_and_corpora", "Neurocognitive_disorder",  "Pandemic_social_impact", "Language_acquisition", "Reading_writing",  "Cancer_illness_studies", "Neuropsychological_assessment", "National_culture", "socio-demographic _research", "Consumer_research", "Aging research", "Psychoanalysis", "social_representation", "Memory_word_retrieval", "Brain_science", "Stroke_cogntivie_function", "Creativity_intelligence", "Community_detection", "concept_biological")

ggraph(stm_corrs, layout = 'fr') +
  geom_edge_link(
    aes(edge_width = weight),
    label_colour = '#fc8d62',
    edge_colour = '#377eb8') +
  geom_node_point(size = 4, colour = 'black')  +
  geom_node_label(
    aes(label = name, size = props),
    colour = 'black',  repel = TRUE, alpha = 0.85) +
  scale_size(range = c(2, 8), labels = scales::percent) +
  labs(size = 'Topic Proportion',  edge_width = 'Topic Correlation') +
  scale_edge_width(range = c(1, 3)) +
  theme_graph()
