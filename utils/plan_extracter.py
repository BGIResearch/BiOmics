def plan_exetract(plans):
    res = []
    for plan in plans:
        res.append(plan.get("type"))
    return res

