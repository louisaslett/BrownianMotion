sim.condbbintersection_ <- function(bm, s_idx, q, t_idx, label) {
  bb.s_idx <- match(bm$t[s_idx], bm$bb.local$t)
  bb.t_idx <- match(bm$t[t_idx], bm$bb.local$t)

  Bs <- bm$W_t[s_idx]
  Bt <- bm$W_tm[t_idx]
  Ws <- bm$bb.local$W_t[bb.s_idx]
  Wt <- bm$bb.local$W_t[bb.t_idx]
  s <- bm$t[s_idx]
  t <- bm$t[t_idx]

  sim.condintersection_(bm$bb.local, bb.s_idx, q, bb.t_idx, label)
  Wq <- bm$bb.local$W_t[bb.s_idx+1]
  Bq <- Bs + (Wq-Ws) + (q-s)/(t-s)*((Bt-Bs)-(Wt-Ws))

  bm$t <- c(bm$t[1:s_idx],
            q,
            bm$t[t_idx:length(bm$t)])
  bm$W_t <- c(bm$W_t[1:s_idx],
              Bq,
              bm$W_t[t_idx:length(bm$W_t)])
  bm$W_tm <- c(bm$W_tm[1:s_idx],
               Bq,
               bm$W_tm[t_idx:length(bm$W_tm)])

  cur.layer <- which(bm$layers$t.l == s & bm$layers$t.u == t)
  i <- nrow(bm$bb.local$layers)-1
  bm$layers <- add_row(bm$layers,
                       type = "intersection-bb",
                       t.l = bm$bb.local$layers$t.l[i],
                       t.u = bm$bb.local$layers$t.u[i],
                       # Ld = Bs + bm$bb.local$layers$Ld[i] - Ws,
                       # Uu = Bs + bm$bb.local$layers$Uu[i] - Ws + ((Bq-Bs)-(Wq-Ws)),
                       # Lu = Bs + bm$bb.local$layers$Ld[i] - Ws + ((Bq-Bs)-(Wq-Ws)),
                       # Ud = Bs + bm$bb.local$layers$Uu[i] - Ws,
                       Ld = Bs - Ws,
                       Uu = Bs - Ws + ((Bq-Bs)-(Wq-Ws)),
                       Lu = NA,
                       Ud = NA,
                       Lu.hard = bm$bb.local$layers$Lu.hard[i], #ifelse(head(c(0, x.new), -1) > x.new, TRUE, FALSE)
                       Ud.hard = bm$bb.local$layers$Ud.hard[i])
  bm$layers <- add_row(bm$layers,
                       type = "intersection-bb",
                       t.l = bm$bb.local$layers$t.l[i+1],
                       t.u = bm$bb.local$layers$t.u[i+1],
                       # Ld = Bq + bm$bb.local$layers$Ld[i+1] - Wq,
                       # Uu = Bq + bm$bb.local$layers$Uu[i+1] - Wq + ((Bt-Bq)-(Wt-Wq)),
                       # Lu = Bq + bm$bb.local$layers$Ld[i+1] - Wq + ((Bt-Bq)-(Wt-Wq)),
                       # Ud = Bq + bm$bb.local$layers$Uu[i+1] - Wq,
                       Ld = Bq - Wq,
                       Uu = Bq - Wq + ((Bt-Bq)-(Wt-Wq)),
                       Lu = NA,
                       Ud = NA,
                       Lu.hard = bm$bb.local$layers$Lu.hard[i+1], #ifelse(head(c(0, x.new), -1) > x.new, TRUE, FALSE)
                       Ud.hard = bm$bb.local$layers$Ud.hard[i+1])
  bm$layers <- bm$layers[-cur.layer,]

  # Sort the aux path since sorting usually happens in the wrapper which was
  # never called from the aux path's perspective ... our main path will be sorted
  # externally
  bm$bb.local$layers <- bm$bb.local$layers[order(bm$bb.local$layers$t.l),]

  add.labels_(bm, "user", q)
  add.labels_(bm, label, q)

  bm
}
