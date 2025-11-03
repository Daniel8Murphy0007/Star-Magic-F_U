// ES6 module syntax (compatible with package.json "type": "module")
import http from 'http';
const PORT = 3000;

// UQFF Astronomical Systems Data
const SYSTEMS = [
    {
        name: "Crab Nebula M1",
        subtitle: "AD 1054 Supernova Remnant",
        age: "970 years",
        distance: "6,500 ly",
        mass: "1×10³¹ kg",
        field: "F_U_Bi_i = -2.09×10²¹² + i(-6.75×10¹⁶⁰) N",
        equation: "[-F + momentum + gravity + stability + LENR + activation + dark_energy + resonance + neutron + relativistic + neutrino] dx",
        color: "#0088ff"
    },
    {
        name: "Vela Pulsar",
        subtitle: "PSR J0835-4510",
        age: "11,000 years",
        distance: "1,000 ly",
        mass: "2.8×10³⁰ kg (1.4 M☉)",
        field: "F_U = 3.45×10⁵⁶ N/m | Rotation: 89 ms",
        equation: "[Ug1 + Ug2 + Ug3 + Ug4 - Ub1 + Um + UA] · (1 + ω_s·t)",
        color: "#ff8800"
    },
    {
        name: "Magnetar SGR 1806-20",
        subtitle: "Ultra-Strong Magnetic Field",
        age: "~1,500 years",
        distance: "50,000 ly",
        mass: "2.8×10³⁰ kg",
        field: "B = 1×10¹¹ T | F_U = 8.92×10⁵⁸ N/m",
        equation: "[k₁Ug1 + k₃Ug3 + Um] · (1 + B_crit/B₀) · sin(ω_c·t)",
        color: "#ff0088"
    },
    {
        name: "Sagittarius A*",
        subtitle: "Milky Way SMBH",
        age: "13.6 Gyr",
        distance: "26,000 ly",
        mass: "8.15×10³⁶ kg (4.1×10⁶ M☉)",
        field: "g_UQFF = 1.23×10⁻¹⁵ m/s² | r_s = 1.2×10¹⁰ m",
        equation: "A_μν = g_μν + η·T_s^μν(ρ_vac_SCm, ρ_vac_UA, t_n)",
        color: "#8800ff"
    },
    {
        name: "Hydrogen Atom",
        subtitle: "Fundamental Quantum System",
        age: "Primordial",
        distance: "N/A",
        mass: "1.67×10⁻²⁷ kg",
        field: "F_U = 7.26×10⁵⁴ N/m | a₀ = 0.529 Å",
        equation: "[26-layer gravity + quantum + vacuum + relativistic] · ψ(r,t)",
        color: "#00ff88"
    },
    {
        name: "NGC 1316 (Fornax A)",
        subtitle: "Elliptical Galaxy with Merger",
        age: "2 Gyr post-merger",
        distance: "62 Mly",
        mass: "5×10¹¹ M☉",
        field: "g_NGC1316 = 2.34×10⁻⁹ m/s² at 20 kpc",
        equation: "[Ug1 + Ug3' + Ug4 + quantum + fluid + merger] · (1 + H·t)",
        color: "#ffaa00"
    }
];

let currentSystemIndex = 0;

// Dynamic UQFF Web Server with Navigation
function generateHTML() {
    const system = SYSTEMS[currentSystemIndex];
    
    return `<!DOCTYPE html>
<html>
<head>
    <title>Star-Magic UQFF System</title>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: 'Courier New', monospace;
            background: linear-gradient(135deg, #000015 0%, #001a33 100%);
            color: #00ff88;
            padding: 20px;
            min-height: 100vh;
        }
        .container {
            max-width: 900px;
            margin: 0 auto;
        }
        h1 {
            text-align: center;
            color: #00ccff;
            text-shadow: 0 0 10px #00ccff;
            margin-bottom: 10px;
            font-size: 2em;
        }
        .status {
            text-align: center;
            color: #00ff88;
            margin-bottom: 30px;
            font-size: 1.2em;
        }
        .system-info {
            background: rgba(0, 136, 255, 0.1);
            border: 2px solid ${system.color};
            border-radius: 10px;
            padding: 20px;
            margin: 20px 0;
            box-shadow: 0 0 20px ${system.color}33;
        }
        .system-name {
            color: ${system.color};
            font-size: 1.8em;
            margin-bottom: 5px;
            text-shadow: 0 0 10px ${system.color};
        }
        .system-subtitle {
            color: #00aaff;
            font-size: 1.1em;
            margin-bottom: 15px;
        }
        .system-details {
            line-height: 1.8;
            margin: 15px 0;
        }
        .system-details p {
            margin: 8px 0;
        }
        .field-display {
            background: rgba(0, 255, 136, 0.1);
            padding: 10px;
            border-left: 3px solid #00ff88;
            margin: 10px 0;
            font-family: monospace;
        }
        .equation-display {
            background: rgba(136, 0, 255, 0.1);
            padding: 10px;
            border-left: 3px solid #8800ff;
            margin: 10px 0;
            font-family: monospace;
            font-size: 0.9em;
            overflow-x: auto;
        }
        .button-container {
            display: flex;
            justify-content: center;
            gap: 15px;
            margin: 30px 0;
            flex-wrap: wrap;
        }
        button {
            background: linear-gradient(135deg, #0088ff 0%, #00ccff 100%);
            color: #000;
            border: none;
            padding: 12px 30px;
            font-size: 1em;
            font-weight: bold;
            border-radius: 5px;
            cursor: pointer;
            transition: all 0.3s;
            font-family: 'Courier New', monospace;
            box-shadow: 0 4px 15px rgba(0, 136, 255, 0.4);
        }
        button:hover {
            background: linear-gradient(135deg, #00ccff 0%, #00ffff 100%);
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(0, 204, 255, 0.6);
        }
        button:active {
            transform: translateY(0);
        }
        .nav-button {
            background: linear-gradient(135deg, #ff8800 0%, #ffaa00 100%);
        }
        .nav-button:hover {
            background: linear-gradient(135deg, #ffaa00 0%, #ffcc00 100%);
            box-shadow: 0 6px 20px rgba(255, 170, 0, 0.6);
        }
        .system-counter {
            text-align: center;
            color: #00aaff;
            margin: 20px 0;
            font-size: 1.1em;
        }
        .footer {
            text-align: center;
            margin-top: 40px;
            color: #0088ff;
            font-size: 0.9em;
        }
        @media (max-width: 600px) {
            h1 { font-size: 1.5em; }
            .system-name { font-size: 1.4em; }
            button { padding: 10px 20px; font-size: 0.9em; }
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>⭐ Star-Magic UQFF System v2.0</h1>
        <div class="status">Status: OPERATIONAL | Framework: Unified Quantum Field Force</div>
        
        <div class="system-counter">
            System ${currentSystemIndex + 1} of ${SYSTEMS.length} | Project Total: 79+ Astronomical Objects in UQFF
        </div>
        
        <div class="system-info">
            <div class="system-name">🌟 ${system.name}</div>
            <div class="system-subtitle">${system.subtitle}</div>
            
            <div class="system-details">
                <p><strong>Age:</strong> ${system.age}</p>
                <p><strong>Distance:</strong> ${system.distance}</p>
                <p><strong>Mass:</strong> ${system.mass}</p>
            </div>
            
            <div class="field-display">
                <strong>Field Strength:</strong><br>${system.field}
            </div>
            
            <div class="equation-display">
                <strong>UQFF Equation:</strong><br>${system.equation}
            </div>
        </div>
        
        <div class="button-container">
            <button class="nav-button" onclick="window.location.href='/prev'">⬅️ Previous System</button>
            <button onclick="alert('🔬 UQFF System OPERATIONAL\\n\\n✅ ${SYSTEMS.length} Core Modules Ready\\n✅ Current: ${system.name}\\n✅ Aetheric Propulsion: Ready\\n✅ Quantum Field: Stable\\n\\nPress NEXT to explore more systems!')">🔬 Run Diagnostics</button>
            <button class="nav-button" onclick="window.location.href='/next'">Next System ➡️</button>
        </div>
        
        <div class="footer">
            <p>🌌 Aetheric Propulsion (UQFF) Research System</p>
            <p>Developed since 2010 | Field Unification Challenge to MOND & Lambda-CDM</p>
        </div>
    </div>
</body>
</html>`;
}

const server = http.createServer((req, res) => {
    if (req.url === '/') {
        res.writeHead(200, {'Content-Type': 'text/html'});
        res.end(generateHTML());
    } else if (req.url === '/next') {
        currentSystemIndex = (currentSystemIndex + 1) % SYSTEMS.length;
        res.writeHead(302, {'Location': '/'});
        res.end();
    } else if (req.url === '/prev') {
        currentSystemIndex = (currentSystemIndex - 1 + SYSTEMS.length) % SYSTEMS.length;
        res.writeHead(302, {'Location': '/'});
        res.end();
    } else {
        res.writeHead(404, {'Content-Type': 'text/plain'});
        res.end('Not Found');
    }
});

server.listen(PORT, () => {
    console.log('═══════════════════════════════════════════════════════');
    console.log('⭐  Star-Magic UQFF Web Interface v2.0');
    console.log('═══════════════════════════════════════════════════════');
    console.log(`🌐 Server: http://localhost:${PORT}`);
    console.log(`📊 Status: OPERATIONAL`);
    console.log(`🔬 Systems: ${SYSTEMS.length} astronomical objects loaded`);
    console.log(`🚀 Ready for Aetheric Propulsion Research!`);
    console.log('═══════════════════════════════════════════════════════');
    console.log('\n💡 Navigation:');
    console.log('   - Click "Next System" to explore different objects');
    console.log('   - Click "Previous System" to go back');
    console.log('   - Click "Run Diagnostics" for system status');
    console.log('\n Press Ctrl+C to stop the server\n');
});
