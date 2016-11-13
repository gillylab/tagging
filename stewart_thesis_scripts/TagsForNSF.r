dir1 = ('/Users/jstewart/Dropbox/Thesis/Dg_Tagging/_CA_Tags/')
setwd(paste(dir1, sep=''))

file = '64004_09/64004-Series.csv'

f = read.csv(file.path(dir1, file), header=T, skip = 3); head(f)

x = f$Depth

w=hist(x, breaks = seq(0, 1400, 50), plot=FALSE);
b = t(t((as.numeric(w$breaks))))
c = t(t((as.numeric(w$counts))))
w2 = cbind(b,c[1:dim(b)[1]])

# need to break it up by day and night

a = strsplit(file, '/', fixed=FALSE)
a.1 = unlist(a)[1]

write.csv(w2, paste(dir1, a.1, '_histo50m.csv', sep=''), row.names=FALSE)