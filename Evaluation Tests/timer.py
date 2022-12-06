import time


class Timer():
    def __init__(self):
        self.start_time = None

    def start(self):
        assert self.start_time == None
        self.start_time = time.perf_counter()

    def stop(self):
        assert self.start_time != None
        tot_time = time.perf_counter() - self.start_time
        self.start_time = None
        return (f"{tot_time:0.3f}")

if __name__ == "__main__":
    t = Timer()
    for _ in range(5):
        t.start()
        time.sleep(2)
        print(t.stop())


