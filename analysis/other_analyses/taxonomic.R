source("analysis/InitialSetup.R")


OTU.tax %>% group_by(Phylum, Class, Order, Family) %>% count(Family) %>%
  ggplot(aes(x = Order, y = n, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  coord_flip()

filter(OTU.tax, OTU == "Otu000001")
OTUs.ds <- OTUsREL[which(design$order > 2),]
OTUs.hw <- OTUsREL[which(design$order <= 2),]
rowSums(OTUs.ds > 0 *1)

species.by.site <- as.data.frame(t(OTUs))
species.by.site <- rownames_to_column(species.by.site, var = "OTU")
otu.table <- full_join(OTU.tax, species.by.site)

otu.table <- otu.table %>% group_by(Phylum) %>% gather(key = sample, value = abundance, -OTU, -Domain, -Phylum, -Class, -Order, -Family, -Genus)
otu.table <- otu.table %>% full_join(rownames_to_column(design, var = "sample")) 
str(otu.table)

otu.table %>% group_by(order) %>% tally(abundance, sort = T) %>% 
  left_join(otu.table) %>% mutate(relabund = abundance/n) %>%
  filter(Phylum != "unclassified") %>%
  ggplot(aes(x = order, y = relabund, fill = Phylum)) + 
  geom_bar(stat = "identity") + 
  coord_flip()
