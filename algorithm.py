def alg(MI_list, gamma, indices, cols):
    t = 0
    while t<1:
        indices.append([])
        cols.append([])
        for x in MI_list[-1][:-1]:
            indices[-1].append(x)
            cols[-1].append(cols[indices.index(x)])
            cols.pop(indices.index(x))
            indices.remove(x)

        t += 1

    return indices, cols
