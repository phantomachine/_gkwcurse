class gkwmod:
    def __init__(self, siggma, betta, thetta, alfa, A, qbar):
        self.siggma = siggma
        self.betta = betta
        self.thetta = thetta  
        self.alfa = alfa
        self.A  = A
        self.qbar = qbar
        description = "Marginal utility, cost, marginal cost functions"
        author = "Timothy Kam (tcy.kam@gmail.com)"

    def utildiff(self, q):
        if self.thetta == 1:
            uprime = 1/(q + self.qbar)
        elif self.thetta > 1 or (self.thetta < 1 and self.thetta > 0):
            uprime = (q + self.qbar)**(-self.thetta)
        else:
            print('Error: THETTA must be positive!')

        return uprime

    def cost(self, q):
        return self.A * q**self.alfa

    def costdiff(self, q):
        return self.alfa * self.A * q**(self.alfa - 1)
