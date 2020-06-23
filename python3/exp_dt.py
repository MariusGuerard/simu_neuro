dt_list = [0.03, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]

for dt in dt_list:
    print(dt)
    exec(compile(open("init.py").read(), "init.py", 'exec'))
