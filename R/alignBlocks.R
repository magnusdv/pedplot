# Family-block alignment built directly from pedigree relations

.alignBlocks = function(input) {
  ids = input$id

  # Remove hidden founders and their parent references before alignment
  hidden = which(input$hidden)
  blockInput = input[!input$hidden, , drop = FALSE]
  blockInput$fid[blockInput$fid %in% ids[hidden]] = 0
  blockInput$mid[blockInput$mid %in% ids[hidden]] = 0

  # Start with earliest possible generations
  ped = .blockPed(blockInput)
  earlyPed = .blockSchedule(ped)
  earlyRows = .blockRows(earlyPed)

  # Also try later generations
  if(nrow(ped$fam)) {
    latePed = .blockSchedule(.blockLate(ped))
    lateRows = .blockRows(latePed)
    nEarly = sum(vapply(earlyRows, function(r)
      sum(vapply(r, function(b) length(b$idx), integer(1L))), integer(1L)))
    nLate = sum(vapply(lateRows, function(r)
      sum(vapply(r, function(b) length(b$idx), integer(1L))), integer(1L)))
    # Keep the earliest improved schedule unless the other removes copies
    if(nLate < nEarly) {
      ped = latePed
      rows = lateRows
    }
    else {
      ped = earlyPed
      rows = earlyRows
    }
  }
  else {
    ped = earlyPed
    rows = earlyRows
  }

  # Order rows and place copies
  rows = .blockPlan(rows, ped)

  # Convert the chosen blocks
  out = .blockOutput(rows, ped, blockInput)

  # Restore hidden-parent details in mating table
  matings = out$matings
  if(length(hidden)) {
    first = blockInput$idx[vapply(ped$kids, `[`, integer(1L), 1L)]
    matings$fid = as.character(input$fid[first])
    matings$mid = as.character(input$mid[first])
    matings$faIdx = match(matings$fid, ids, nomatch = 0L)
    matings$moIdx = match(matings$mid, ids, nomatch = 0L)
    matings$faHidden = matings$faIdx %in% hidden
    matings$moHidden = matings$moIdx %in% hidden
  }

  # Finalise alignment data
  nodes = out$nodes
  dup = .blockDupEdges(nodes)
  edges = if(length(dup$edge)) rbind(out$edges, dup) else out$edges
  edges$edge = seq_len(nrow(edges))
  row.names(edges) = NULL

  align = list(
    input = input, nodes = nodes, matings = matings, edges = edges,
    meta = list(hidden = ids[input$hidden],
                duplicated = unique.default(ids[dup$child]),
                maxlev = max(nodes$y)))

  class(align) = "pedalign"
  align
}


# Move matings between safe rows only when this joins person occurrences
.blockSchedule = function(ped) {
  fam = ped$fam
  if(!nrow(fam))
    return(ped)

  # Store every family in which each person is a parent
  parent = c(fam$faIdx[fam$faIdx > 0L], fam$moIdx[fam$moIdx > 0L])
  mating = c(fam$mating[fam$faIdx > 0L], fam$mating[fam$moIdx > 0L])
  parentMats = split.default(mating, factor(parent, levels = seq_len(ped$n)))

  # Skip the passes when every person already occupies one generation
  mixed = FALSE
  for(i in seq_len(ped$n)) {
    lev = fam$gen[parentMats[[i]]]
    birth = ped$familyOf[i]
    if(birth > 0L)
      lev = c(lev, fam$gen[birth] + 1L)
    if(length(lev) > 1L && any(lev != lev[1L])) {
      mixed = TRUE
      break
    }
  }
  ceiling = max(ped$depth) - 1L

  # Use two fixed down-and-up passes instead of an open-ended search
  ord = if(mixed) order(fam$gen, fam$mating) else integer()
  laterMats = if(mixed) lapply(ped$kids, function(k)
    unlist(parentMats[k], use.names = FALSE)) else list()

  for(f in rep(c(ord, rev.default(ord)), 2L)) {
    current = fam$gen[f]
    par = c(fam$faIdx[f], fam$moIdx[f])
    par = par[par > 0L]
    kids = ped$kids[[f]]

    # Bound this mating between parental births and child matings
    parBirth = ped$familyOf[par]
    lower = if(any(parBirth > 0L))
      max(1L, fam$gen[parBirth[parBirth > 0L]] + 1L) else 1L
    later = laterMats[[f]]
    upper = if(length(later)) min(fam$gen[later] - 1L) else ceiling

    # Try only rows which can join an affected person's other occurrence
    cand = current
    parRest = vector("list", length(par))
    for(j in seq_along(par)) {
      p = par[j]
      other = parentMats[[p]]
      rest = fam$gen[other[other != f]]
      birth = ped$familyOf[p]
      if(birth > 0L)
        rest = c(rest, fam$gen[birth] + 1L)
      parRest[[j]] = rest
      cand = c(cand, rest)
    }

    kidRest = vector("list", length(kids))
    for(j in seq_along(kids)) {
      rest = fam$gen[parentMats[[kids[j]]]]
      kidRest[[j]] = rest
      cand = c(cand, rest - 1L)
    }
    cand = cand[cand >= lower & cand <= upper]

    # Count new person-generation occurrences for every candidate row
    cost = numeric(length(cand))
    for(rest in parRest) if(length(rest))
      cost = cost + !cand %in% rest
    for(rest in kidRest) if(length(rest))
      cost = cost + !(cand + 1L) %in% rest
    best = which.min((cost * (ceiling + 1L) + abs(cand - current)) *
                     (ceiling + 1L) + cand)

    # Accept only a strict occurrence reduction
    if(cost[best] < cost[1L])
      fam$gen[f] = cand[best]
  }

  # Rebuild primary generations from the chosen family rows
  ped$fam = fam
  founder = ped$familyOf == 0L
  ped$depth[!founder] = fam$gen[ped$familyOf[!founder]] + 1L
  ped$depth[founder] = 1L
  for(i in which(founder)) if(length(parentMats[[i]]))
    ped$depth[i] = min(fam$gen[parentMats[[i]]])

  ped
}

# Move short branches down without changing any ancestry order
.blockLate = function(ped) {
  # Collect all matings in which each person is a parent
  fam = ped$fam
  nFam = nrow(fam)
  parent = c(fam$faIdx[fam$faIdx > 0L], fam$moIdx[fam$moIdx > 0L])
  mating = c(fam$mating[fam$faIdx > 0L], fam$mating[fam$moIdx > 0L])
  parentMats = split.default(mating, factor(parent, levels = seq_len(ped$n)))

  # Put every terminal family on the lowest available row
  latest = rep(max(ped$depth) - 1L, nFam)

  # Process descendants before the family in which each person was born
  for(i in rev.default(ped$order)) {
    birth = ped$familyOf[i]
    later = parentMats[[i]]
    if(birth > 0L && length(later))
      latest[birth] = min(latest[birth], min(latest[later]) - 1L)
  }

  # Derive primary generations from the delayed family rows
  ped$fam$gen = latest
  founder = ped$familyOf == 0L

  # Place nonfounders one row below their birth family
  ped$depth[!founder] = latest[ped$familyOf[!founder]] + 1L

  # Place each used founder at its earliest delayed mating
  for(i in which(founder)) {
    mating = parentMats[[i]]
    if(length(mating))
      ped$depth[i] = min(latest[mating])
  }
  ped
}


# Build parent indexes, safe generations and one row per mating
.blockPed = function(input) {
  # Reject empty input before constructing indexes
  n = nrow(input)
  if(!n)
    stop2("No individuals")

  fa = .parentIdx(input$fid, input$id, "father")
  mo = .parentIdx(input$mid, input$id, "mother")

  # Reject one impossible parent pair explicitly
  if(any(fa > 0L & fa == mo))
    stop2("Same father and mother: ", input$id[fa > 0L & fa == mo][1L])

  # Assign generations from parent-child links only
  dep = .blockDepth(fa, mo)
  depth = dep$depth
  child = which(fa > 0L | mo > 0L)

  if(length(child)) {
    # Keep children with the same ordered parent pair together
    key = fa[child] * (n + 1L) + mo[child]
    group = match(key, unique.default(key))
    nFam = max(group)

    # Preserve the first child and input sibling order
    first = match(seq_len(nFam), group)
    kids = unname(split.default(child, factor(group, levels = seq_len(nFam))))
    familyOf = integer(n)
    familyOf[child] = group

    # Read each parent pair from its first child
    faFam = fa[child[first]]
    moFam = mo[child[first]]
    faGen = moGen = integer(nFam)
    faGen[faFam > 0L] = depth[faFam[faFam > 0L]]
    moGen[moFam > 0L] = depth[moFam[moFam > 0L]]

    # Start every mating at its deeper parent's row
    fam = list2DF(list(mating = seq_len(nFam), faIdx = faFam,
                       moIdx = moFam, gen = pmax(faGen, moGen)))
  }
  else {
    # Preserve the same table shapes for pedigrees without matings
    familyOf = integer(n)
    kids = list()
    fam = list2DF(list(mating = integer(), faIdx = integer(),
                       moIdx = integer(), gen = integer()))
  }

  # Give every person a stable structural ordering value
  anchor = as.numeric(seq_len(n))
  for(i in dep$order) {
    # Inherit the mean parental anchor when parents are known
    par = c(fa[i], mo[i])
    par = par[par > 0L]
    if(length(par))
      anchor[i] = mean.default(anchor[par])
  }

  list(n = n, fa = fa, mo = mo, depth = depth, order = dep$order,
       familyOf = familyOf, fam = fam, kids = kids, anchor = anchor)
}


# Find longest ancestry depths in one topological pass
.blockDepth = function(fa, mo) {
  # Count the unresolved parents of every person
  n = length(fa)
  indeg = (fa > 0L) + (mo > 0L)
  child = c(which(fa > 0L), which(mo > 0L))
  parent = c(fa[fa > 0L], mo[mo > 0L])

  # Store each child under both known parents
  kids = split.default(child, factor(parent, levels = seq_len(n)))

  # Start with all founders in input order
  queue = integer(n)
  roots = which(indeg == 0L)
  queue[seq_along(roots)] = roots
  head = 1L
  tail = length(roots)
  depth = rep(1L, n)
  topo = integer(n)
  done = 0L

  # Visit each parent-child link exactly once
  while(head <= tail) {
    # Remove the next person from the fixed queue
    p = queue[head]
    head = head + 1L
    done = done + 1L
    topo[done] = p

    # Release children after their last parent is visited
    for(ch in kids[[p]]) {
      depth[ch] = max(depth[ch], depth[p] + 1L)
      indeg[ch] = indeg[ch] - 1L
      if(indeg[ch] == 0L) {
        # Append newly available children without growing the queue
        tail = tail + 1L
        queue[tail] = ch
      }
    }
  }

  if(done < n)
    stop2("Impossible pedigree: parent-child cycle")

  # Return both depths and the reusable topological order
  list(depth = depth, order = topo)
}


# Build linear spouse blocks and add every remaining person once
.blockRows = function(ped) {
  # Split people and matings once instead of scanning each generation
  fam = ped$fam
  nGen = max(ped$depth)
  rows = vector("list", nGen)
  people = split.default(seq_len(ped$n), factor(ped$depth, levels = seq_len(nGen)))
  two = which(fam$faIdx > 0L & fam$moIdx > 0L)
  one = which(xor(fam$faIdx > 0L, fam$moIdx > 0L))
  twoByGen = split.default(two, factor(fam$gen[two], levels = seq_len(nGen)))
  oneByGen = split.default(one, factor(fam$gen[one], levels = seq_len(nGen)))

  for(g in seq_along(rows)) {
    # Merge all two-parent family blocks on this row into trails
    two = twoByGen[[g]]
    blocks = .blockTrails(fam$faIdx[two], fam$moIdx[two], fam$mating[two])

    # Orient complete trails without changing their mating links
    if(length(blocks))
      blocks = lapply(blocks, .blockOrient, fam = fam, anchor = ped$anchor)

    # Find people already supplied by spouse trails
    used = unique.default(unlist(lapply(blocks, `[[`, "idx"), use.names = FALSE))
    one = oneByGen[[g]]
    onePar = c(fam$faIdx[one], fam$moIdx[one])
    onePar = onePar[onePar > 0L]

    # Add primary people and parents from one-parent families
    need = unique.default(c(people[[g]], onePar))
    alone = need[!need %in% used]

    # Give every singleton its own ordering group
    nGroup = length(blocks)
    if(length(alone)) {
      alone = alone[order(ped$anchor[alone], alone)]

      # Create singleton blocks in one pre-sized list
      singles = lapply(seq_along(alone), function(i)
        list(idx = alone[i], mating = integer(), group = nGroup + i))
      blocks = c(blocks, singles)
    }

    rows[[g]] = blocks
  }

  rows
}


# Return flat occurrence indexes for one ordered row
.blockInfo = function(row) {
  len = vapply(row, function(b) length(b$idx), integer(1L))
  n = sum(len)
  primary = unlist(lapply(row, function(b) b$primary), use.names = FALSE)

  list(idx = unlist(lapply(row, `[[`, "idx"), use.names = FALSE),
       block = rep(seq_along(row), len), local = sequence.default(len),
       primary = primary, x = seq_len(n) - (n + 1L) / 2)
}


# Order rows from their parent families and add copies to prevent crossings
.blockPlan = function(rows, ped) {
  fam = ped$fam
  nFam = nrow(fam)
  matPos = rep(NA_real_, nFam)
  personPos = rep(NA_real_, ped$n)
  one = which(xor(fam$faIdx > 0L, fam$moIdx > 0L))
  oneByGen = split.default(one, factor(fam$gen[one], levels = seq_along(rows)))
  people = split.default(seq_len(ped$n), factor(ped$depth, levels = seq_along(rows)))

  for(g in seq_along(rows)) {
    ownPeople = people[[g]]
    # Start with no selected occurrence on this row
    rr = lapply(rows[[g]], function(b) {
      b$primary = rep(FALSE, length(b$idx))
      b
    })
    nBlock = length(rr)

    if(nBlock > 1L) {
      # Locate each block from birth families, then copies and anchors
      key = rep(NA_real_, nBlock)
      anch = first = numeric(nBlock)
      for(b in seq_len(nBlock)) {
        ids = rr[[b]]$idx
        birth = ped$familyOf[ids]
        own = ped$depth[ids] == g & birth > 0L
        target = matPos[birth[own]]
        target = target[is.finite(target)]
        if(length(target))
          key[b] = stats::median.default(target)
        else {
          copy = personPos[ids[ped$depth[ids] < g]]
          copy = copy[is.finite(copy)]
          if(length(copy))
            key[b] = stats::median.default(copy)
        }
        anch[b] = mean.default(ped$anchor[unique.default(ids)])
        first[b] = min(ids)
      }

      # Put unanchored blocks on the same scale as parent positions
      known = is.finite(key)
      if(any(known)) {
        arank = rank(anch, ties.method = "first")
        lim = range.default(key[known])
        fallback = sum(lim)/2 +
          (arank - (nBlock + 1L)/2) * max(lim[2L] - lim[1L], 1) / max(nBlock - 1L, 1L)
        key[!known] = fallback[!known]
      }
      else
        key = rank(anch, ties.method = "first")

      # Reorder complete blocks with stable structural ties
      ord = order(key, anch, first, seq_len(nBlock))
      rr = rr[ord]
    }

    # Select one own-generation occurrence for every person
    info = .blockInfo(rr)
    for(i in ownPeople) {
      cand = which(info$idx == i)
      birth = ped$familyOf[i]
      if(birth > 0L && is.finite(matPos[birth]))
        cand = cand[which.min(abs(info$x[cand] - matPos[birth]))]
      else
        cand = cand[1L]
      rr[[info$block[cand]]]$primary[info$local[cand]] = TRUE
    }

    # Turn primary children towards their ordered birth families
    rr = .blockFace(rr, ped, matPos, personPos)

    # Rank birth families by the mating directly above
    birthFam = unique.default(ped$familyOf[ownPeople])
    birthFam = birthFam[birthFam > 0L]
    famRank = integer(nFam)
    if(length(birthFam)) {
      birthFam = birthFam[order(matPos[birthFam], birthFam)]
      famRank[birthFam] = seq_along(birthFam)
    }

    # Repair rigid-boundary inversions with a bounded insertion pass
    nOwn = length(ownPeople)
    for(repair in seq_len(nOwn * nOwn)) {
      info = .blockInfo(rr)
      take = which(info$primary & ped$familyOf[info$idx] > 0L)
      label = famRank[ped$familyOf[info$idx[take]]]
      bad = which(diff(label) < 0L)
      if(!length(bad))
        break
      bad = bad[1L] + 1L
      person = info$idx[take[bad]]
      higher = take[seq_len(bad - 1L)][label[seq_len(bad - 1L)] > label[bad]]
      insert = min(info$block[higher])
      current = info$block[take[bad]]

      # Retire all occurrences of this person
      occ = which(info$idx == person)
      for(k in occ)
        rr[[info$block[k]]]$primary[info$local[k]] = FALSE

      # Reuse an existing earlier copy before adding another
      cand = occ[info$block[occ] < insert]
      if(length(cand)) {
        chosen = cand[length(cand)]
        rr[[info$block[chosen]]]$primary[info$local[chosen]] = TRUE
      }
      else if(length(rr[[current]]$idx) == 1L && !length(rr[[current]]$mating)) {
        # Move an unconstrained singleton instead of duplicating it
        rr[[current]]$primary = TRUE
        moving = rr[current]
        rr = rr[-current]
        if(current < insert)
          insert = insert - 1L
        rr = append(rr, moving, after = insert - 1L)
      }
      else {
        # Keep the mating block and add one birth-side copy
        group = max(vapply(rr, `[[`, integer(1L), "group")) + 1L
        copy = list(idx = person, mating = integer(), group = group,
                    primary = TRUE)
        rr = append(rr, list(copy), after = insert - 1L)
      }
    }

    # Save compact row positions for ordering the next generation
    rows[[g]] = rr
    info = .blockInfo(rr)
    own = which(info$primary)
    personPos[info$idx[own]] = info$x[own]
    for(b in seq_along(rr)) if(length(rr[[b]]$mating)) {
      pos = info$x[info$block == b]
      matPos[rr[[b]]$mating] = (pos[-length(pos)] + pos[-1L]) / 2
    }

    # Locate one-parent matings at their plotted parent occurrence
    for(m in oneByGen[[g]]) {
      par = max(fam$faIdx[m], fam$moIdx[m])
      cand = which(info$idx == par)
      primary = cand[info$primary[cand]]
      chosen = if(length(primary)) primary[1L] else cand[1L]
      matPos[m] = info$x[chosen]
    }
  }

  rows
}


# Orient rigid blocks from primary birth targets before cosmetic ties
.blockFace = function(row, ped, matPos, personPos) {
  info = .blockInfo(row)
  for(b in seq_along(row)) {
    ids = row[[b]]$idx
    loc = which(row[[b]]$primary & ped$familyOf[ids] > 0L)
    target = matPos[ped$familyOf[ids[loc]]]

    # Prefer monotone birth targets before placing copies outwards
    inv1 = if(length(target) > 1L) sum(pmax(0, -diff(target))) else 0
    inv2 = if(length(target) > 1L) sum(pmax(0, diff(target))) else 0
    pos = info$x[info$block == b]
    old = personPos[ids]
    dist1 = sum(abs(pos - old), na.rm = TRUE)
    dist2 = sum(abs(rev.default(pos) - old), na.rm = TRUE)
    mat = row[[b]]$mating
    fa = ped$fam$faIdx[mat]
    left1 = sum(ids[-length(ids)] == fa)
    left2 = sum(ids[-1L] == fa)
    reverse = if(inv1 != inv2) inv2 < inv1 else
      if(dist1 != dist2) dist2 < dist1 else left2 > left1
    if(reverse) for(nm in c("idx", "mating", "primary"))
      row[[b]][[nm]] = rev.default(row[[b]][[nm]])
  }
  row
}


# Choose trail direction from structure, spouse side and input order
.blockOrient = function(block, fam, anchor) {
  # Singletons have no direction to choose
  idx = block$idx
  if(length(idx) < 2L)
    return(block)

  # Prefer the direction agreeing best with structural order
  pos = seq_along(idx)
  score = rank(anchor[idx], ties.method = "average")
  cost1 = sum(abs(score - pos))
  cost2 = sum(abs(rev.default(score) - pos))

  # Reverse for a better structural direction, then input order
  ri = rev.default(idx)
  fst = which(idx != ri)[1L]
  reverse = cost2 < cost1 || (cost1 == cost2 && !is.na(fst) && ri[fst] < idx[fst])

  if(reverse) {
    # Reverse vertices and matching edges together
    block$idx = ri
    block$mating = rev.default(block$mating)
  }
  block
}


# Decompose each spouse component into the fewest open trails
.blockTrails = function(u, v, mating) {
  # An empty row has no spouse components
  if(!length(u))
    return(list())

  # Use compact vertex numbers within this generation
  person = sort.int(unique.default(c(u, v)))
  u = match(u, person)
  v = match(v, person)

  # Build one adjacency list for component discovery
  edge = seq_along(u)
  n = length(person)
  adj = split.default(c(edge, edge), factor(c(u, v), levels = seq_len(n)))
  comp = integer(n)
  nComp = 0L
  queue = integer(n)

  # Label spouse-connected components once
  for(root in seq_len(n)) {
    if(comp[root] > 0L)
      next
    nComp = nComp + 1L

    # Start a fixed-size breadth-first queue
    head = tail = 1L
    queue[1L] = root
    comp[root] = nComp
    while(head <= tail) {
      # Visit every incident mating once from each endpoint
      a = queue[head]
      head = head + 1L
      for(e in adj[[a]]) {
        b = if(u[e] == a) v[e] else u[e]
        if(comp[b] == 0L) {
          # Label and queue each newly reached spouse
          tail = tail + 1L
          queue[tail] = b
          comp[b] = nComp
        }
      }
    }
  }

  groups = split.default(edge, factor(comp[u], levels = seq_len(nComp)))

  # Count the minimum number of trails before allocating output
  nTrail = vapply(groups, function(ii) {
    ends = c(u[ii], v[ii])
    degree = tabulate(match(ends, ends))
    nOdd = sum(degree %% 2L == 1L)
    if(nOdd) nOdd %/% 2L else 1L
  }, integer(1L))
  ans = vector("list", sum(nTrail))
  out = 0L

  for(cc in seq_along(groups)) {
    # Compact the vertices again within this component
    ii = groups[[cc]]
    vert = sort.int(unique.default(c(u[ii], v[ii])))
    uu = match(u[ii], vert)
    vv = match(v[ii], vert)
    degree = tabulate(c(uu, vv), nbins = length(vert))

    # Odd vertices become the endpoints of open trails
    odd = which(degree %% 2L == 1L)

    if(length(odd)) {
      # A dummy vertex turns all minimum trails into one circuit
      dummy = length(vert) + 1L
      uu = c(uu, rep(dummy, length(odd)))
      vv = c(vv, odd)
      tag = c(mating[ii], integer(length(odd)))

      # Traverse every real and dummy edge exactly once
      z = .blockEuler(uu, vv, dummy)
      keep = z$vertex != dummy
      runs = split.default(which(keep), cumsum(!keep)[keep])

      # Removing dummy links recovers the open spouse trails
      for(pp in runs) {
        # Map compact vertices back to pedigree indexes
        out = out + 1L
        ans[[out]] = list(idx = person[vert[z$vertex[pp]]],
                          mating = tag[z$edge[pp[-1L]]], group = cc)
      }
    }
    else {
      # An even component is one closed trail with one repeated person
      z = .blockEuler(uu, vv, 1L)

      # Retain the repeated start needed to open the circuit
      out = out + 1L
      ans[[out]] = list(idx = person[vert[z$vertex]],
                        mating = mating[ii][z$edge[-1L]], group = cc)
    }
  }

  ans
}


# Return a deterministic Euler circuit and its incoming edges
.blockEuler = function(u, v, start) {
  # Store each undirected edge under both endpoints
  m = length(u)
  n = max(c(u, v))
  edge = seq_len(m)
  adj = split.default(c(edge, edge), factor(c(u, v), levels = seq_len(n)))
  nextEdge = rep(1L, n)
  used = logical(m)

  # Hierholzer's stack uses fixed storage of one slot per edge
  stackV = stackE = integer(m + 2L)
  stackV[1L] = start
  top = 1L
  circuitV = circuitE = integer(m + 1L)
  out = 0L

  while(top > 0L) {
    # Advance past edges consumed on another visit
    at = stackV[top]
    inc = adj[[at]]
    p = nextEdge[at]
    while(p <= length(inc) && used[inc[p]])
      p = p + 1L
    nextEdge[at] = p

    if(p <= length(inc)) {
      # Follow the next unused edge and remember its arrival
      e = inc[p]
      used[e] = TRUE
      nextEdge[at] = p + 1L
      top = top + 1L
      stackV[top] = if(u[e] == at) v[e] else u[e]
      stackE[top] = e
    }
    else {
      # Add a dead end while unwinding the current trail
      out = out + 1L
      circuitV[out] = at
      circuitE[out] = stackE[top]
      top = top - 1L
    }
  }

  # Reverse the unwound circuit into traversal order
  take = rev.default(seq_len(out))
  list(vertex = circuitV[take], edge = circuitE[take])
}


# Fit monotone block shifts and interpolate unconstrained blocks
.blockFit = function(target, n) {
  known = which(is.finite(target))
  if(!length(known))
    return(numeric(n))
  fitted = stats::isoreg(known, target[known])$yf
  if(length(known) == 1L)
    return(rep(fitted, n))
  stats::approx(known, fitted, seq_len(n), rule = 2)$y
}

.blockRefine = function(x, rows, ped, faOcc, moOcc, primaryOcc, snap = 0.2) {
  nFam = nrow(ped$fam)
  if(!nFam)
    return(x)

  len = unlist(lapply(rows, function(r)
    vapply(r, function(b) length(b$idx), integer(1L))), use.names = FALSE)
  nBlock = length(len)
  block = rep(seq_len(nBlock), len)
  first = cumsum(c(1L, len[-nBlock]))
  last = first + len - 1L

  # Minimise child distances from their parental centre
  D = diag(1e-8, nBlock)
  d = numeric(nBlock)
  centre = matrix(0, nFam, nBlock)
  centreB = numeric(nFam)

  for(m in seq_len(nFam)) {
    par = c(faOcc[m], moOcc[m])
    par = par[par > 0L]
    ch = primaryOcc[ped$kids[[m]]]
    w = 1/length(ch)
    px = mean.default(x[par])

    # Store the exact family-centre equation for snapping
    u = unique.default(c(block[par], block[ch]))
    centre[m, u] =
      tabulate(match(block[par], u), length(u))/length(par) -
      tabulate(match(block[ch], u), length(u))/length(ch)
    centreB[m] = mean.default(x[ch]) - px

    for(k in ch) {
      u = unique.default(c(block[par], block[k]))
      a = tabulate(match(block[par], u), length(u))/length(par)
      a[match(block[k], u)] = a[match(block[k], u)] - 1
      r = px - x[k]

      D[u, u] = D[u, u, drop = FALSE] + w * tcrossprod(a)
      d[u] = d[u] - w * a * r
    }
  }

  # Preserve block order and minimum spacing within each row
  row = rep(seq_along(rows), lengths(rows))
  adj = which(diff(row) == 0L)
  cmat = matrix(0, length(adj), nBlock)
  b = numeric(length(adj))

  if(length(adj)) {
    j = seq_along(adj)
    cmat[cbind(j, adj)] = -1
    cmat[cbind(j, adj + 1L)] = 1
    b = 1 - (x[first[adj + 1L]] - x[last[adj]])
    shift = quadprog::solve.QP(D, d, t(cmat), b)$solution
  }
  else
    shift = solve(D, d)

  # Make already-near vertical family stems exactly vertical
  hit = abs(drop(centre %*% shift) - centreB) <= snap
  if(any(hit))
    shift = quadprog::solve.QP(
      D, d, t(rbind(centre[hit, , drop = FALSE], cmat)),
      c(centreB[hit], b), meq = sum(hit))$solution

  x + shift[block]
}

# Convert ordered blocks directly to pedalign tables
.blockOutput = function(rows, ped, input) {
  # Count occurrences before allocating node columns
  fam = ped$fam
  nFam = nrow(fam)
  nOcc = sum(vapply(rows, function(rr)
    sum(vapply(rr, function(b) length(b$idx), integer(1L))), integer(1L)))

  idx = gen = integer(nOcc)
  x = numeric(nOcc)
  isPrimary = logical(nOcc)
  primaryOcc = integer(ped$n)
  faOcc = moOcc = integer(nFam)
  one = which(xor(fam$faIdx > 0L, fam$moIdx > 0L))
  oneByGen = split.default(one, factor(fam$gen[one], levels = seq_along(rows)))
  matByGen = split.default(seq_len(nFam), factor(fam$gen, levels = seq_along(rows)))
  at = 0L

  # Place each row from parent centres while preserving rigid blocks
  for(g in seq_along(rows)) {
    rr = rows[[g]]
    len = vapply(rr, function(b) length(b$idx), integer(1L))
    nBlock = length(rr)
    base = cumsum(c(0L, len[-nBlock]))
    birthQ = copyQ = rep(NA_real_, nBlock)

    # Calculate the preferred start of every constrained block
    for(b in seq_len(nBlock)) {
      block = rr[[b]]
      ids = block$idx
      loc = which(block$primary & ped$familyOf[ids] > 0L)
      if(length(loc)) {
        target = numeric(length(loc))
        for(j in seq_along(loc)) {
          m = ped$familyOf[ids[loc[j]]]
          par = c(faOcc[m], moOcc[m])
          par = par[par > 0L]
          target[j] = mean.default(x[par])
        }
        birthQ[b] = stats::median.default(target - (loc - 1L)) - base[b]
      }

      # Use earlier primary occurrences only when no child constrains the block
      old = which(ped$depth[ids] < g & primaryOcc[ids] > 0L)
      if(length(old))
        copyQ[b] = stats::median.default(x[primaryOcc[ids[old]]] - (old - 1L)) - base[b]
    }

    # Use copy targets for blocks without a birth target
    birthQ[!is.finite(birthQ)] = copyQ[!is.finite(birthQ)]
    q = .blockFit(birthQ, nBlock)
    start = base + q
    rowFirst = at + 1L

    for(b in seq_len(nBlock)) {
      block = rr[[b]]
      # Fill one complete block without growing any vector
      len = length(block$idx)
      ii = at + seq_len(len)
      idx[ii] = block$idx
      gen[ii] = g
      x[ii] = start[b] + seq_len(len) - 1L
      isPrimary[ii] = block$primary
      primaryOcc[block$idx[block$primary]] = ii[block$primary]

      # Record the exact parent copies used by each mating
      if(length(block$mating)) {
        left = ii[-length(ii)]
        right = ii[-1L]
        mat = block$mating

        # Recover father and mother nodes from actual parent indexes
        leftFather = idx[left] == fam$faIdx[mat]
        faOcc[mat] = left + !leftFather
        moOcc[mat] = right - !leftFather
      }
      at = at + len
    }

    # Locate one-parent matings before placing their children
    row = rowFirst:at
    for(m in oneByGen[[g]]) {
      par = max(fam$faIdx[m], fam$moIdx[m])
      cand = row[idx[row] == par]
      primary = cand[isPrimary[cand]]
      chosen = if(length(primary)) primary[1L] else cand[1L]
      if(fam$faIdx[m] > 0L) faOcc[m] = chosen else moOcc[m] = chosen
    }
  }

  x = .blockRefine(x, rows, ped, faOcc, moOcc, primaryOcc)

  # Centre the completed plot without changing relative positions
  x = x - (max(x) + min(x))/2
  byPerson = split.default(seq_len(nOcc), factor(idx, levels = seq_len(ped$n)))

  # Allocate explicit copy numbers for every occurrence
  copy = integer(nOcc)

  # Number the selected primary first and all other copies by position
  for(i in seq_len(ped$n)) {
    # Sort non-primary copies from top to bottom and left to right
    jj = byPerson[[i]]
    other = jj[jj != primaryOcc[i]]
    if(length(other) > 1L)
      other = other[order(gen[other], x[other], other)]

    copy[primaryOcc[i]] = 1L
    copy[other] = seq_along(other) + 1L
  }

  nodes = list2DF(list(
    node = seq_len(nOcc), idx = input$idx[idx], id = input$id[idx], copy = copy,
    x = x, y = length(rows) - gen + 1L,  gen = gen, sex = input$sex[idx],
    primary = isPrimary))

  # Build mating rows with the same fields as the existing aligners
  faHidden = moHidden = logical(nFam)
  hasFa = fam$faIdx > 0L
  hasMo = fam$moIdx > 0L

  # Read hidden status only for known parent indexes
  faHidden[hasFa] = input$hidden[fam$faIdx[hasFa]]
  moHidden[hasMo] = input$hidden[fam$moIdx[hasMo]]
  fid = mid = rep("0", nFam)

  # Restore external parent labels for the public tables
  fid[hasFa] = input$id[fam$faIdx[hasFa]]
  mid[hasMo] = input$id[fam$moIdx[hasMo]]
  single = xor(!hasFa | faHidden, !hasMo | moHidden)

  # Keep node references even when a founder is removed later
  matings = list2DF(list(
    mating = fam$mating, fid = fid, mid = mid, faIdx = fam$faIdx,
    moIdx = fam$moIdx, faNode = faOcc, moNode = moOcc, faHidden = faHidden,
    moHidden = moHidden, singleParent = single, gen = fam$gen))

  # Add every child relation and every two-parent mating once
  child = unlist(ped$kids, use.names = FALSE)
  childMat = rep(fam$mating, lengths(ped$kids))

  # Select every two-parent family for one spouse edge
  nChild = length(child)
  spouseMat = fam$mating[hasFa & hasMo]

  # Mark matings whose parents share an ancestor
  par = unique.default(c(fam$faIdx[spouseMat], fam$moIdx[spouseMat]))
  anc = vector("list", ped$n)
  for(i in par)
    anc[[i]] = c(i, .blockAncestor(i, ped$mo, ped$fa))

  double = vapply(spouseMat, function(m)
    any(anc[[fam$faIdx[m]]] %in% anc[[fam$moIdx[m]]]), logical(1L))

  edges = list2DF(list(
    edge = seq_len(nChild + length(spouseMat)),
    type = c(rep("child", nChild), rep("spouse", length(spouseMat))),
    from = c(childMat, faOcc[spouseMat]),
    to = c(primaryOcc[child], moOcc[spouseMat]),
    mating = c(childMat, rep(NA_integer_, length(spouseMat))),
    child = c(input$idx[child], rep(NA_integer_, length(spouseMat))),
    code = c(rep("", nChild), ifelse(double, "2", "1"))))

  # Duplicate edges and final table cleanup are shared downstream
  list(nodes = nodes, matings = matings, edges = edges)
}

# Return all ancestors using a fixed queue
.blockAncestor = function(me, mom, dad) {
  seen = logical(length(mom))
  queue = integer(length(mom))
  queue[1L] = me
  head = tail = 1L

  while(head <= tail) {
    par = c(mom[queue[head]], dad[queue[head]])
    head = head + 1L
    par = par[par > 0L]
    par = par[!seen[par]]

    if(length(par)) {
      seen[par] = TRUE
      queue[tail + seq_along(par)] = par
      tail = tail + length(par)
    }
  }

  which(seen)
}

# Link adjacent copies after ordering duplicated nodes once
.blockDupEdges = function(nodes) {
  idx = nodes$idx
  node = nodes$node
  ord = which(tabulate(idx)[idx] > 1L)
  ord = ord[order(idx[ord], nodes$x[ord], nodes$gen[ord], node[ord])]

  left = ord[-length(ord)]
  right = ord[-1L]
  same = idx[left] == idx[right]
  left = left[same]
  right = right[same]
  n = length(left)

  list2DF(list(edge = rep(NA_integer_, n),
               type = rep("duplicate", n),
               from = node[left],
               to = node[right],
               mating = rep(NA_integer_, n),
               child = idx[right],
               code = rep("", n)))
}
