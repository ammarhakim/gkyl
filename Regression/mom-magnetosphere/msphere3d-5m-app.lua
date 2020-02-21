local buildApp = require("App.Apps.MSPHERE")

local app = buildApp {
   tEnd=2700,
}

app:run()
