class BimObject:
    def __init__(self, bim_path):
        self._bim_path = open(bim_path, "r")


    def dosomething(self):

        print(self._bim_path.readline())

            # print(f.readline())
            # f.seek(96)
            # print(f.readline())
            #
            # cumulative_seek = 0
            # indexer = {}
            # for reader_line in f:
            #     line = f.readline()
            #
            #     if len(line) > 0:
            #         _, rsid, _, _, _, _ = line.split()
            #         indexer[rsid] = cumulative_seek
            #         cumulative_seek += len(line)
            #
            #     # if len(line) > 0:
            #     #     _, rsid, _, _, _, _ = line.split()
            #     #     print(rsid)
            #     #     print(len(line))
            #
            # print(indexer)
