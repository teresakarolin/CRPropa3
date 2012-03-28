from mpc import *
from pylab import *
import sys
import ROOT

gamma = logspace(6, 14, 200)

def get_rates(interaction, id, gamma, count=10000):
    D = {}
    c = Candidate()
    c.current.setId(id)
    c.current.setLorentzFactor(gamma)

    for i in range(count):
        c.clearInteractionStates()
        interaction.process(c)
        state = InteractionState()
        c.getInteractionState("mpc::PhotoDisintegration", state)
        D.setdefault(state.channel, []).append(state.distance / Mpc)
        
    for channel in D.keys():
        l = D[channel]
        hist = ROOT.TH1F('', '', 50, 0, max(l))
        for element in l:
            hist.Fill(element)
        f = ROOT.TF1('f1', 'expo')
        hist.Fit(f, "q")
        l = -f.GetParameter(1)
        p = hist.GetEntries() / float(count)
        exclu = l * p # exclusive decay constant for this channel
        D[channel] = max(0, exclu) # take 0 if slope is negative
    
    return D

def get_rates_per_channel(interaction, id):
    D = {}
    for iE in range(200):
        if iE % 20 != 0:
            continue
        D2 = get_rates(interaction, id, gamma[iE])
        for channel in D2.keys():
            D.setdefault(channel, zeros(200))[iE] = D2[channel]
    return D

def get_data_rates_per_channel(table, id):
    rates = {}
    for i in xrange(len(table)):
        if int(table[i][0]) != id:
            continue
        rates[table[i][1]] = table[i][2:]
    return rates

def plot_channel(rates_simulated, rates_data, id, channel):
    figure()
    if channel in rates_data:
        plot(gamma, rates_data[channel], 'r', label="Data")
    if channel in rates_simulated:
        plot(gamma, rates_simulated[channel], 'k+', label="Simulated")
    legend(loc='lower right')
    text(0.1, 0.85, 'Nucleus ' + parse_id(id) + '\nDecay Channel ' + parse_channel(channel), transform=gca().transAxes)
    xlabel('Lorentzfactor $\gamma$')
    ylabel('Rate [1/Mpc]')
    loglog()
    ylim(1e-10, 1e2)
    grid()
    savefig('PhotoDisintegration_' + str(id) + '_' + str(channel) + '.png', bbox_inches='tight')

def get_id(a, z):
    return 1e9 + a * 1e6 + z * 1e3

def parse_id(id):
    z = ((id - 1000000000) % 1000000) // 1000
    a = (id - 1000000000) // 1000000
    return 'Z=%i, A=%i' % (z, a)

def parse_channel(c):
    s = '%06d' % c
    d = list(map(int, s))
    s = 'n, ' * d[0] + 'p, ' * d[1] + 'H$^2$, ' * d[2] + 'H$^3$, ' * d[3] + 'He$^3$, ' * d[4] + 'He$^4$, ' * d[5]
    return s[0:-2]

interaction = PhotoDisintegration()
table = genfromtxt(getDataPath('/PhotoDisintegration/pd_table.txt'))
id = 1004002000

if len(sys.argv) >= 3:
    a = int(sys.argv[1])
    z = int(sys.argv[2])
    id = get_id(a, z)

print 'Plotting Disintegration Rates for', id
rates_simulated = get_rates_per_channel(interaction, id)
rates_data = get_data_rates_per_channel(table, id)
channels = set.union(set(rates_simulated.keys()), set(rates_data.keys()))
for channel in channels:
    plot_channel(rates_simulated, rates_data, id, channel)

# plot channel multiplicity
print 'Plotting Channel Multiplicity'
multi = zeros((27, 31))
for z in range(1, 27):
    for n in range(1, 31):
        multi[z][n] = len(get_data_rates_per_channel(table, get_id(z + n, z)))

fig = figure()
ax = fig.add_subplot(111)
mmulti = ma.masked_array(multi, multi==0)
im = ax.imshow(mmulti, aspect='equal', interpolation='nearest', origin='lower')
cbar = fig.colorbar(im, orientation='horizontal')
cbar.set_label('Photo-Disintegration Channels')
ax.set_xlabel('Neutrons')
ax.set_ylabel('Protons')
ax.grid()
fig.savefig('PhotoDisintegration_multiplicity.png',bbox_inches='tight')


