// Ready-made experiments: two molecules, settings, what to do and what comes out.
// Every entry was run in tests/experiments-screen.mjs (two runs each, the pull done by the
// same spring as the mouse) and gave the stated products both times; `energy` is the change
// of the potential energy seen there, in eV.
//   products: what the line "now in the scene" reads when it has worked (charges are not shown there)
//   pull: what to take hold of and where to lead it; without it the molecules react on their own
//   etemp: electronic temperature the experiment needs (K), default 300
//   unpaired: number of unpaired electrons, if not the scene's own (see scene.js)
//   speed: fs of simulation per second
export const EXPERIMENTS = [
  { key: 'neutral', products: '2 H₂O', group: 'Acids and bases', title: 'Acid + base → water', first: 'h3o+', second: 'oh-', etemp: 1500, speed: 200,
    todo: 'Press “Start” and wait: the ions attract each other on their own.', expect: 'A proton passes from H₃O⁺ to OH⁻, giving two water molecules. This is neutralization; the energy drops by about 6.6 eV.' },
  { key: 'acid-nh3', products: 'H₂O + NH₄', group: 'Acids and bases', title: 'Acid + ammonia → ammonium', first: 'h3o+', second: 'nh3', speed: 200,
    todo: 'Press “Start” and wait.', expect: 'Ammonia takes the proton: H₂O and the ammonium ion NH₄⁺. The energy drops by about 1.9 eV.' },
  { key: 'nh4-oh', products: 'NH₃ + H₂O', group: 'Acids and bases', title: 'Ammonium + base → ammonia', first: 'nh4+', second: 'oh-', etemp: 1500, speed: 200,
    todo: 'Press “Start” and wait.', expect: 'Hydroxide takes a proton from ammonium: NH₃ and water. This is why an ammonium salt smells of ammonia once you add a base.' },
  { key: 'naoh-hf', products: 'NaF·H₂O', group: 'Acids and bases', title: 'Caustic soda + acid → salt and water', first: 'naoh', second: 'hf', speed: 60,
    pull: 'Press “Start”, grab the hydrogen of the HF molecule and bring it to the oxygen of NaOH.', expect: 'The proton moves to the oxygen: sodium fluoride and water, which stay stuck together. The energy drops by about 2.4 eV.' },
  { key: 'hf-oh', products: 'F + H₂O', group: 'Acids and bases', title: 'HF + hydroxide → water and fluoride', first: 'hf', second: 'oh-', speed: 60,
    pull: 'Press “Start”, grab the oxygen of OH⁻ and bring it to the hydrogen of HF.', expect: 'Water and the F⁻ ion. The energy drops by about 3.3 eV.' },
  { key: 'cyanide', products: 'H₂O + HCN', group: 'Acids and bases', title: 'Cyanide + acid → hydrogen cyanide', first: 'h3o+', second: 'cn-', etemp: 1500, speed: 200,
    todo: 'Press “Start” and wait.', expect: 'Cyanide takes the proton: HCN and water. This is why cyanides must never be acidified.' },

  { key: 'soda', products: 'NaHCO₃', group: 'Carbon dioxide', title: 'Base + CO₂ → baking soda', first: 'naoh', second: 'co2', speed: 60,
    pull: 'Press “Start”, grab the carbon of CO₂ and bring it to the oxygen of NaOH.', expect: 'The oxygen of the base attaches to the carbon and CO₂ bends: sodium bicarbonate NaHCO₃, baking soda. The energy drops by 1.6–1.9 eV.' },
  { key: 'hco3', products: 'HCO₃', group: 'Carbon dioxide', title: 'Hydroxide + CO₂ → bicarbonate', first: 'oh-', second: 'co2', speed: 60,
    pull: 'Press “Start”, grab the carbon of CO₂ and bring it to the oxygen of OH⁻.', expect: 'The HCO₃⁻ ion. This is how alkaline solutions absorb carbon dioxide. The energy drops by about 2.7 eV.' },

  { key: 'lih-hf', products: 'LiF + H₂', group: 'Hydrogen', title: 'Lithium hydride + acid → hydrogen', first: 'lih', second: 'hf', speed: 60,
    pull: 'Press “Start”, grab the hydrogen of the HF molecule and bring it to the hydrogen of LiH.', expect: 'Two hydrogens, one with an excess of electrons and the other with a shortage, join into H₂ and fly off; LiF is left behind. The energy drops by about 1.6 eV.' },
  { key: 'bh4-acid', products: 'BH₃ + H₂ + H₂O', group: 'Hydrogen', title: 'Borohydride + acid → hydrogen', first: 'bh4-', second: 'h3o+', etemp: 1500, speed: 200,
    todo: 'Press “Start” and wait.', expect: 'The proton of the acid meets a hydridic hydrogen of the borohydride: H₂, borane BH₃ and water.' },
  { key: 'li-f2', products: 'Li₂F₂', group: 'Hydrogen', title: 'Lithium burns in fluorine', first: 'li2', second: 'f2', etemp: 1500, speed: 200,
    todo: 'Press “Start” and wait.', expect: 'The electrons of lithium jump to fluorine, the molecules fall apart and reassemble into the ionic square Li₂F₂. The energy drops by about 10 eV: the most violent reaction in the list.' },

  { key: 'li-o2', products: 'Li₂O₂', group: 'Oxygen', title: 'Lithium is oxidized by oxygen', first: 'li2', second: 'o2', speed: 200,
    todo: 'Press “Start” and wait.', expect: 'The electrons of lithium go over to oxygen, giving the ionic rhombus Li₂O₂, lithium peroxide. The energy drops by about 2.8 eV. The two unpaired electrons of oxygen do not go anywhere: spin is conserved.' },
  { key: 'h2-o2', products: 'H₂ + O₂', group: 'Oxygen', title: 'Oxyhydrogen without a spark', first: 'h2', second: 'o2', speed: 300,
    todo: 'Press “Start” and wait. Try pushing the molecules together by hand.', expect: 'Nothing. Hydrogen and oxygen do not react at room temperature: for burning to start you need fragments of molecules, which a spark provides. In six test runs of 4 ps each the molecules only bounced off each other.' },

  { key: 'hh', products: 'H₂', group: 'Radicals and spin', title: 'Two hydrogen atoms → H₂', first: 'h', second: 'h', speed: 100,
    todo: 'Press “Start” and wait.', expect: 'Each atom has one electron, and the spins are opposite. The atoms attract and form an H₂ molecule; the energy drops by about 5.2 eV. This is the chemical bond in its purest form.' },
  { key: 'hh-parallel', products: '2 H', unpaired: 2, group: 'Radicals and spin', title: 'The same atoms, spins parallel', first: 'h', second: 'h', speed: 100,
    todo: 'Press “Start” and try to bring the atoms together by hand. Then, under “Settings”, change “Unpaired e⁻” from 2 to 0.', expect: 'No bond: two electrons with the same spin cannot occupy one orbital, so the atoms only repel. Switch the number of unpaired electrons to 0 and they join at once.' },
  { key: 'li-f', products: 'LiF', group: 'Radicals and spin', title: 'Lithium atom + fluorine atom → salt', first: 'li', second: 'f', speed: 100,
    todo: 'Press “Start” and wait.', expect: 'The single valence electron of lithium jumps to fluorine and the ions attract: LiF. The energy drops by about 4.1 eV.' },
  { key: 'ch3-h', products: 'CH₄', group: 'Radicals and spin', title: 'Methyl + hydrogen atom → methane', first: 'ch3', second: 'h', speed: 60,
    pull: 'Press “Start”, grab the hydrogen atom and bring it to the carbon.', expect: 'The two unpaired electrons form a C–H bond and the flat methyl becomes a tetrahedron: methane. The energy drops by about 5.7 eV.' },
  { key: 'ch3-ch3', products: 'C₂H₆', group: 'Radicals and spin', title: 'Two methyls → ethane', first: 'ch3', second: 'ch3', speed: 60,
    pull: 'Press “Start”, grab the carbon of one methyl and bring it to the carbon of the other.', expect: 'A C–C bond: ethane. This is how chain reactions terminate, when two radicals meet. The energy drops by about 5.4 eV.' },
  { key: 'h-f2', products: 'HF + F', group: 'Radicals and spin', title: 'A hydrogen atom breaks a fluorine molecule', first: 'f2', second: 'h', speed: 60,
    pull: 'Press “Start”, grab the hydrogen atom and bring it to either fluorine.', expect: 'HF forms, and the second fluorine flies off as an atom with an unpaired electron, ready to break the next molecule. This is one link of a chain reaction; the energy drops by about 3.6 eV.' },

  { key: 'adduct', products: 'BH₃NH₃', group: 'Boron', title: 'Borane + ammonia → adduct', first: 'bh3', second: 'nh3', speed: 60,
    pull: 'Press “Start”, grab the nitrogen and bring it to the boron.', expect: 'The lone pair of nitrogen enters the empty orbital of boron, and the borane bends out of plane. The energy drops by about 0.8 eV.' },
  { key: 'bf3-nh3', products: 'BF₃NH₃', group: 'Boron', title: 'BF₃ + ammonia → adduct', first: 'bf3', second: 'nh3', speed: 60,
    pull: 'Press “Start”, grab the nitrogen and bring it to the boron.', expect: 'The same as with borane: a B–N bond, the energy drops by about 0.8 eV.' },
  { key: 'diborane', products: 'B₂H₆', group: 'Boron', title: 'Two boranes → diborane', first: 'bh3', second: 'bh3', speed: 60,
    pull: 'Press “Start”, grab the boron of one molecule and bring it to the boron of the other.', expect: 'Two hydrogens become bridges between the boron atoms: B₂H₆. The energy drops by about 0.4 eV.' },
  { key: 'bh3-co', products: 'BH₃CO', group: 'Boron', title: 'Borane + carbon monoxide', first: 'bh3', second: 'co', speed: 60,
    pull: 'Press “Start”, grab the carbon of CO and bring it to the boron.', expect: 'CO attaches to boron through its carbon: borane carbonyl H₃B–CO. The energy drops by about 1.75 eV.' },
  { key: 'diborane-heat', products: '2 BH₃', group: 'Boron', title: 'Heat diborane', first: 'b2h6', second: '', temperature: 2000, speed: 300,
    todo: 'Press “Start” and wait: the thermostat is set to 2000 K.', expect: 'The bridges break and diborane splits into two boranes.' },

  { key: 'lif-lif', products: 'Li₂F₂', group: 'Salts', title: 'Two LiF molecules → seed of a crystal', first: 'lif', second: 'lif', speed: 60,
    pull: 'Press “Start”, grab the fluorine of one molecule and bring it to the lithium of the other.', expect: 'The ions line up in the square Li₂F₂, plus to minus: this is how a salt crystal begins. The energy drops by about 1.3 eV.' },
  { key: 'naf-water', products: 'NaF·H₂O', group: 'Salts', title: 'Salt + water', first: 'naf', second: 'h2o', speed: 60,
    pull: 'Press “Start”, grab the oxygen of water and bring it to the sodium.', expect: 'Water sticks to the sodium ion through its oxygen: the first step of dissolving. The energy drops by 0.5–0.6 eV; a single water molecule does not break the salt apart.' },
];

// Things one might expect to react but which do not here, and why (shown nowhere yet).
export const NOT_REACTING = ['h2o + co2', 'nh3 + co2', 'h2 + f2', 'h2 + o3 at room temperature'];
